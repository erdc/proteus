import copy
from math import *
from EGeometry import *
from LinearAlgebra import *
from LinearSolvers import *
from MeshTools import *
from FemTools import *
from QuadTools import *
from HamiltonJacobiTimeIntegrationTools import *
from NonlinearSolvers import *
import femIntegrals
"""
A module for discretizing Hamilton-Jacobi equations.
"""

## \defgroup HamiltonJacobi HamiltonJacobi
#
# A module for discretizing Hamilton-Jacobi equations.
#
# @{

class HamiltonJacobiCoefficients:
    def __init__(self,
                 mass='nonlinear',
                 hamiltonian='nonlinear'):
        self.mass=mass
        self.hamiltonian
    def evaluate(self,
                 t,
                 x,
                 u,grad_u,
                 m,dm,
                 h,dh,
                 rh):
        pass

class OneLevelHamiltonJacobi(NonlinearEquation):
    
    """
    A class for  finite element discretizations of  Hamilton-Jacobi
    equations on a single spatial mesh.

    Objects of this type take the initial-boundary value problem for

    u_t + H(u,\grad u) = 0

    and turn it into the discrete (nonlinear or linear) algebraic
    problem

    F(U) = 0

    where F and U have dimension self.dim. The Jacobian of F or an
    approximation for  it may also be  provided.

    The NonlinearEquation interface is

    self.dim
    getResidual(u,r)
    getJacobian(jacobian)

    The rest of the functions in this class are either private functions
    or return various other pieces of information.
    """
    
    """keys for each integral in the residual"""
    integralKeys = ['m','h','stab','numDiff']
    
    """keys for each integral in the residual (element boundaries) --these are for scalar transport"""
    elementBoundaryIntegralKeys = ['f']
    
    """keys for each coefficient of the pde"""
    coefficientKeys = ['m','dm','h','dh']
    
    def __init__(self,
                 u,
                 testSpace,
                 matType,
                 dofBoundaryConditions,
                 coefficients,
                 quadrature,
                 elementBoundaryQuadrature,
                 stabilization='1',
                 shockCapturing='1',
                 shockCapturingDiffusion=0.1,
                 TimeIntegrationClass=BackwardEuler,
                 calculateElementBoundaryValues = False,
		 tOrder=1):

        """
        Allocate storage and initialize some variables.

        u -- a ScalarFiniteElementFunction object
        testSpace -- a FiniteElementSpace object,
                     dim(testSpace) = dim(u.femSpace)
        dofBoundaryConditions -- a DOFBoundaryConditions object
        coefficients -- a HamiltonJacobiCoefficients object
        quadrature -- a dictionary of quadrature rules for each integral
        stabilization -- an optional string
        shockCapturing -- an optional string

        The constructure sets the input arguments, calculates some
        dimensions, and allocates some storage Dimension subscripts
        are:

        _global -- per physical domain
        _element -- per element
        
        Storage is divided into quantities at:
        quadrature points - nElements_global x
                            nQuadraturePoints_element x
                            dim(quantity)
        element nodes - nElements_global x nDOF_element x dim(quantity)
        global nodes - nDOF_global
        """
	#keys for each integral in the residual"""
	self.integralKeys = ['m','h','stab','numDiff']
	
	#keys for each integral in the residual (element boundaries) --these are for scalar transport"""
	self.elementBoundaryIntegralKeys = ['f']
	
	#keys for each coefficient of the pde"""
	self.coefficientKeys = ['m','dm','h','dh']
	
        #
        #set the objects describing the method and boundary conditions
        #
        self.u = u
        self.trialSpace = self.u.femSpace
        self.matType = matType
        self.mesh = self.u.femSpace.elementMaps.mesh
        self.testSpace = testSpace
        self.dirichletConditions = dofBoundaryConditions
        self.coefficients = coefficients
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.shockCapturingDiffusion = shockCapturingDiffusion
        self.numericalFlux = None
        self.testUpdateIntegrals=False#True#False
        self.testUpdateJacobianIntegrals=False#True
        self.testCalculateShape=False#True
        self.stabilizationIsNonlinear = not (self.coefficients.hamiltonian == 'linear')
        self.dirichletNodeSetList=None
        #determine if we need element boundary storage
        self.elementBoundaryIntegrals = calculateElementBoundaryValues
        #
        #calculate some dimensions
        #
        self.nSpace_global = u.femSpace.referenceFiniteElement.referenceElement.dim
        #we'll use the max dof per element in anticipation of p-refinement
        #this should probably be calculated by the FiniteElementSpace object
        self.nDOF_element = self.u.femSpace.max_nDOF_element
        self.nFreeDOF_global = self.dirichletConditions.nFreeDOF_global
        #
        NonlinearEquation.__init__(self,self.nFreeDOF_global)
        #
        #get the union of all quadrature points
        #
        quadraturePointSet = set()
        for I in self.integralKeys:
            quadraturePointSet |= set(
                [(p[X],p[Y],p[Z]) for p in quadrature[I].points])
        self.nQuadraturePoints_element = len(quadraturePointSet)
        self.nQuadraturePoints_global = (self.mesh.nElements_global*
                                         self.nQuadraturePoints_element)
        #
        #build a dictionary at each quadrature point which
        #contains a dictionary of weights for each integral
        #
        # e.g. quadratureWeightDict[p]['u'] is the weight at p for the
        # u integral
        #
        #initialize all weights to 0
        quadratureWeightDict={}
        for k,p in enumerate(quadraturePointSet):
            quadratureWeightDict[p]={}
            for I in self.integralKeys:
                quadratureWeightDict[p][I]=0.0
        #set the nonzero weights
        for I in self.integralKeys:
            for w,p in zip(quadrature[I].weights,quadrature[I].points):
                quadratureWeightDict[(p[X],p[Y],p[Z])][I]=w
        #create point and weight arrays
        self.quadraturePoints = Numeric.zeros((self.nQuadraturePoints_element,3),
                                              Numeric.Float)
        for k,p in enumerate(quadraturePointSet):
            self.quadraturePoints[k][:]=p
        self.quadratureWeights = {}
        for I in self.integralKeys:
            self.quadratureWeights[I] = Numeric.zeros(
                (self.nQuadraturePoints_element,),Numeric.Float)
            for k,p in enumerate(quadraturePointSet):
                self.quadratureWeights[I][k] = quadratureWeightDict[p][I]
        #
        #Repeat the same thing for the element boundary quadrature
        #
        #First calculate the union of all element quadrature points.
        #
        elementBoundaryQuadraturePointSet = set()
        for I in self.elementBoundaryIntegralKeys:
            elementBoundaryQuadraturePointSet |= set(
                [(p[X],p[Y],p[Z]) for p in elementBoundaryQuadrature[I].points])
        self.nElementBoundaryQuadraturePoints_elementBoundary = len(elementBoundaryQuadraturePointSet)
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global*
                                                        self.mesh.nElementBoundaries_element*
                                                        self.nElementBoundaryQuadraturePoints_elementBoundary)
        #
        #Now build a dictionary at each element quadrature point which
        #contains the weights for each integral
        #
        # e.g. quadratureWeightDict[p]['m'] is the weight at p for the
        # mass integral
        #
        #initialize all weights to 0
        elementBoundaryQuadratureWeightDict={}
        for k,p in enumerate(elementBoundaryQuadraturePointSet):
            elementBoundaryQuadratureWeightDict[p]={}
            for I in self.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureWeightDict[p][I]=0.0
        #set the nonzero weights
        for I in self.elementBoundaryIntegralKeys:
            for w,p in zip(elementBoundaryQuadrature[I].weights,elementBoundaryQuadrature[I].points):
                elementBoundaryQuadratureWeightDict[(p[X],p[Y],p[Z])][I]=w
        #
        # Now create the desired point and weight arrays
        #
        self.elementBoundaryQuadraturePoints = Numeric.zeros((self.nElementBoundaryQuadraturePoints_elementBoundary,3),
                                              Numeric.Float)
        for k,p in enumerate(elementBoundaryQuadraturePointSet):
            self.elementBoundaryQuadraturePoints[k][:]=p
        self.elementBoundaryQuadratureWeights = {}
        for I in self.elementBoundaryIntegralKeys:
            self.elementBoundaryQuadratureWeights[I] = Numeric.zeros(
                (self.nElementBoundaryQuadraturePoints_elementBoundary,),Numeric.Float)
            for k,p in enumerate(elementBoundaryQuadraturePointSet):
                self.elementBoundaryQuadratureWeights[I][k] = elementBoundaryQuadratureWeightDict[p][I]
        #
        #element quantities
        #
        self.e={}
        self.e['diameter'] = Numeric.zeros(
            (self.mesh.nElements_global,),Numeric.Float)
        for eN in range(self.mesh.nElements_global):
            self.e['diameter'][eN] = self.mesh.elementList[eN].diameter
        #time stepping 
        self.T=0.0
        #
        #quadrature point quantities nElements_global x nQuadraturePoints_element x dim object
        #
        self.q={}
        self.scalars_quadrature = ['u',
                                   'm','dm','mt','dmt','dx_m',
                                   'dx_h','h','rh',
                                   'tau','dtau','dx_stab',
                                   'numDiff','dnumDiff','dx_numDiff',
                                   'mt*dx_m',
                                   'dmt*dx_m',
                                   'h*dx_h',
                                   'tau*dx_stab',
                                   'numDiff*dx_numDiff',
                                   'pdeResidual',
                                   'det(J)',
                                   'cfl']
        for k in self.scalars_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,self.nQuadraturePoints_element),
            Numeric.Float)
        self.vectors_quadrature = ['grad(u)',
                                   'grad(phi)',
                                   'dh',
                                   'dhdu',
                                   'dh*dx_h']
        for k in self.vectors_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),
            Numeric.Float)
        self.q['x'] = Numeric.zeros((self.mesh.nElements_global,
                                     self.nQuadraturePoints_element,
                                     3),
                                    Numeric.Float)
        self.tensors_quadrature = ['J',
                                   'inverse(J)']
        for k in self.tensors_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global,self.nSpace_global),
            Numeric.Float)
        self.shape_quadrature = ['v','w']
        if self.stabilization != None:
            self.shape_quadrature += ['dpdeResidual','Lstar*w']
        for k in self.shape_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nDOF_element),Numeric.Float)
        self.gradient_shapeGradient_quadrature = ['grad(u)*grad(w)']
        
        for k in self.gradient_shapeGradient_quadrature:
            self.q[k]=Numeric.zeros((self.mesh.nElements_global,
                                       self.nQuadraturePoints_element,
                                       self.nDOF_element,
                                       self.nSpace_global,self.nSpace_global),Numeric.Float)
        self.shapeGradient_quadrature = ['grad(v)','grad(w)']
        for k in self.shapeGradient_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nDOF_element,
             self.nSpace_global),Numeric.Float)
        self.shape_shape_quadrature = ['v*w']
        for k in self.shape_shape_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nDOF_element,self.nDOF_element),Numeric.Float)
        self.shape_shapeGradient_quadrature = ['v*grad(w)','grad(v)*w']
        for k in self.shape_shapeGradient_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nDOF_element,self.nDOF_element,
             self.nSpace_global),Numeric.Float)
        self.shapeGradient_shapeGradient_quadrature = ['grad(v)*grad(w)',]
        for k in self.shapeGradient_shapeGradient_quadrature:
            self.q[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.nQuadraturePoints_element,
                 self.nDOF_element,self.nDOF_element,
                 self.nSpace_global,self.nSpace_global),
                Numeric.Float)
        #
        #element boundary quadrature point quanties (global) -- nElementBoundaries_global x nElementBoundaryQuadraturePoints_elementBoundary x dim object
        #
        self.ebq_global = {}
        self.ebq = {}
        if self.elementBoundaryIntegrals:
            #element boundary quantities
            self.vectors_elementBoundaryQuadrature_global = ['n']
            for k in self.vectors_elementBoundaryQuadrature_global: self.ebq_global[k]=Numeric.zeros(
                (self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),
                Numeric.Float)
            self.ebq_global['x'] = Numeric.zeros(
                (self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),
                Numeric.Float)
            #
            #element boundary quadrature point quantities (unique to each element)-- nElements_global x nElementBoundaries_element x nElementBoundaryQuadraturePoints_elementBoundary x dim object
            #
            self.scalars_elementBoundaryQuadrature = ['u',
                                                      'sqrt(det(g))']
            for k in self.scalars_elementBoundaryQuadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global, self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary),
                Numeric.Float)
            self.vectors_elementBoundaryQuadrature = ['n',
                                                      'grad(u)']
            for k in self.vectors_elementBoundaryQuadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global),
                Numeric.Float)
            self.points_elementBoundaryQuadrature = ['x','hat(x)','bar(x)']
            for k in self.points_elementBoundaryQuadrature:
                self.ebq[k] = Numeric.zeros((self.mesh.nElements_global,
                                               self.mesh.nElementBoundaries_element,
                                               self.nElementBoundaryQuadraturePoints_elementBoundary,
                                               3),
                                              Numeric.Float)
            self.tensors_elementBoundaryQuadrature = ['inverse(J)']
            for k in self.tensors_elementBoundaryQuadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nSpace_global,
                 self.nSpace_global),
                Numeric.Float)
            self.ebq['g'] = Numeric.zeros((self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_element,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary,
                                           max(1,self.nSpace_global-1),
                                           max(1,self.nSpace_global-1)),
                                          Numeric.Float)
            self.shape_ElementBoundaryQuadrature = ['v','w']
            for k in self.shape_quadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_element),Numeric.Float)
            self.gradient_shapeGradient_normal_elementBoundaryQuadrature = ['grad(u)*n']
            for k in self.gradient_shapeGradient_normal_elementBoundaryQuadrature:
                self.ebq[k]=Numeric.zeros((self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_element,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary,
                                           self.nDOF_element,
                                           self.nSpace_global,self.nSpace_global),Numeric.Float)
            self.shapeGradient_elementBoundaryQuadrature = ['grad(v)']
            for k in self.shapeGradient_elementBoundaryQuadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_element,
                 self.nSpace_global),Numeric.Float)
            self.shape_shape_elementBoundaryQuadrature = ['v*w']
            for k in self.shape_shape_elementBoundaryQuadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_element,self.nDOF_element),Numeric.Float)
            self.shape_shapeGradient_elementBoundaryQuadrature = ['grad(v)*w']
            for k in self.shape_shapeGradient_elementBoundaryQuadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.mesh.nElementBoundaries_element,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_element,self.nDOF_element,
                 self.nSpace_global),Numeric.Float)
        #
        # element residual and Jacobian storage
        #
        self.elementResidual = Numeric.zeros(
            (self.mesh.nElements_global,
             self.nDOF_element),
            Numeric.Float)
        if self.testUpdateIntegrals:
            self.restemp = Numeric.array(self.elementResidual)
        self.elementJacobian = Numeric.zeros(
            (self.mesh.nElements_global,
             self.nDOF_element,
             self.nDOF_element),
            Numeric.Float)
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp = Numeric.array(self.elementJacobian)
        #
        # Build the node connectivity lists for solvers
        #
        # TODO: fix this for DG
        #
        self.mesh.buildNodeStarList()
        self.freeNodeStarList=[]
        for n in range(self.mesh.nNodes_global):
            self.freeNodeStarList.append([])
        for n in range(self.mesh.nNodes_global):
            if self.dirichletConditions.global2freeGlobal.has_key(n):
                N = self.dirichletConditions.global2freeGlobal[n]
                self.mesh.nodeStarList[n].sort()
                for m in self.mesh.nodeStarList[n]:
                    if self.dirichletConditions.global2freeGlobal.has_key(m):
                        M = self.dirichletConditions.global2freeGlobal[m]
                        self.freeNodeStarList[N].append(M)
        print "updating local to global mapping matrix"
        self.updateLocal2Global()
        print "building time inegration object"
        print TimeIntegrationClass
        self.timeIntegration = TimeIntegrationClass(self)
        print "Updating Quadrature"
        self.updateQuadrature()
    def setInitialConditions(self,getInitialConditions,T=0.0):
        interpolationValues = Numeric.zeros((self.u.femSpace.elementMaps.mesh.nElements_global,
                                             self.u.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints),
                                            Numeric.Float)
        for eN in range(self.u.femSpace.elementMaps.mesh.nElements_global):
            for k in range(self.u.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                interpolationValues[eN,k] = getInitialConditions.uOfXT(self.u.femSpace.interpolationPoints[eN,k],T)
        self.u.projectFromInterpolationConditions(interpolationValues)
        #Load the Dirichlet conditions 
        for dofN,g in self.dirichletConditions.DOFBoundaryConditionDict.iteritems():
            self.u.dof[dofN] = g(self.dirichletConditions.DOFBoundaryPointDict[dofN],self.T)
        self.updateCoefficients()
    def initializeTimeIntegration(self):
        self.updateCoefficients()
        #mwf added to make sure RK schemes get right initial conditions
        self.timeIntegration.setInitialStageValues()
        self.timeIntegration.updateTimeHistory()
    def getResidual(self,u,r):
        """
        Calculate the element residuals and add in to the global residual
        """
        r.flat[:]=0.0
        #Load the Dirichlet conditions 
        for dofN,g in self.dirichletConditions.DOFBoundaryConditionDict.iteritems():
            self.u.dof[dofN] = g(self.dirichletConditions.DOFBoundaryPointDict[dofN],self.T)
        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        self.updateCoefficients()
        self.updateElementResidual()
        #Sum up contributions to the global residual
#          for eN in range(self.mesh.nElements_global):
#              localDOF = self.l2g[eN]['freeLocal']
#              freeGlobalDOF = self.l2g[eN]['freeGlobal']
#              for i,I in zip(localDOF,freeGlobalDOF):
#                  r[I] += self.elementResidual[eN,i]
#          for eN in range(self.mesh.nElements_global):
#              for ii in range(self.l2g['nFreeDOF'][eN]):
#                  i=self.l2g['freeLocal'][eN,ii]
#                  I=self.l2g['freeGlobal'][eN,ii]
#                  r[I] += self.elementResidual[eN,i]
        femIntegrals.updateGlobalResidualFromElementResidual(self.mesh.nElements_global,
                                                             self.nDOF_element,
                                                             self.l2g['nFreeDOF'],
                                                             self.l2g['freeLocal'],
                                                             self.l2g['freeGlobal'],
                                                             self.elementResidual,
                                                             r);
        return r
    def getJacobian(self,jacobian):
	if self.matType == 'csr':
	    femIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
					  jacobian)
	else:
	    jacobian.flat[:]=0.0
	self.updateElementJacobian()
	#self.updateElementBoundaryJacobian()
	if self.matType == 'csr':
	    self.getJacobian_CSR(jacobian)
	else:
	    self.getJacobian_dense(jacobian)
	return jacobian
    def getJacobian_dense(self,jacobian):
	femIntegrals.updateGlobalJacobianFromElementJacobian_dense(self.mesh.nElements_global,
								   self.nDOF_element,
								   self.nFreeDOF_global,
								   self.l2g['nFreeDOF'],
								   self.l2g['freeLocal'],
								   self.l2g['freeGlobal'],
								   self.elementJacobian,
								   jacobian)
#         if self.numericalFlux != None:
#             femIntegrals.updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(self.mesh.nInteriorElementBoundaries_global,
# 											   self.nDOF_element,
# 											   self.nElementBoundaryQuadraturePoints_elementBoundary,
# 											   self.mesh.nElementBoundaries_element,
# 											   self.nFreeDOF_global,
# 											   self.mesh.interiorElementBoundariesArray,
# 											   self.mesh.elementBoundaryElementsArray,
# 											   self.mesh.elementBoundaryLocalElementBoundariesArray,
# 											   self.l2g['nFreeDOF'],
# 											   self.l2g['freeLocal'],
# 											   self.l2g['freeGlobal'],
# 											   self.fluxJacobian,
# 											   self.ebq['w'],
# 											   jacobian)
#             femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(self.mesh.nExteriorElementBoundaries_global,
# 											   self.nDOF_element,
# 											   self.nElementBoundaryQuadraturePoints_elementBoundary,
# 											   self.mesh.nElementBoundaries_element,
# 											   self.nFreeDOF_global,
# 											   self.mesh.exteriorElementBoundariesArray,
# 											   self.mesh.elementBoundaryElementsArray,
# 											   self.mesh.elementBoundaryLocalElementBoundariesArray,
# 											   self.l2g['nFreeDOF'],
# 											   self.l2g['freeLocal'],
# 											   self.l2g['freeGlobal'],
# 											   self.fluxJacobian,
# 											   self.ebq['w'],
# 											   jacobian)
#         elif self.fluxBoundaryConditions == 'outFlow':
# 	    femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(self.mesh.nExteriorElementBoundaries_global,
# 											   self.nDOF_element,
# 											   self.nElementBoundaryQuadraturePoints_elementBoundary,
# 											   self.mesh.nElementBoundaries_element,
# 											   self.nFreeDOF_global,
# 											   self.mesh.exteriorElementBoundariesArray,
# 											   self.mesh.elementBoundaryElementsArray,
# 											   self.mesh.elementBoundaryLocalElementBoundariesArray,
# 											   self.l2g['nFreeDOF'],
# 											   self.l2g['freeLocal'],
# 											   self.l2g['freeGlobal'],
# 											   self.fluxJacobian,
# 											   self.ebq['w'],
# 											   jacobian)
	return jacobian
    def getJacobian_CSR(self,jacobian):
        """
        Add in the element jacobians to the global jacobian
        """
 #         if self.numericalFlux != None:
#              for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
#                  ebN = self.mesh.interiorElementBoundariesArray[ebNI]
#                  left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
#                  left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
#                  for ii in range(self.l2g['nFreeDOF'][left_eN_global]):
#                      for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                          jacIndex = self.csrRowIndeces[left_eN_global,ii] + self.csrColumnOffsets_eb[ebN,0,0,ii,jj]
#                          femIntegrals.zeroJacobian(jacIndex,
#                                                    jacobian)
#                      for jj in range(self.l2g['nFreeDOF'][right_eN_global]):
#                          jacIndex = self.csrRowIndeces[left_eN_global,ii] + self.csrColumnOffsets_eb[ebN,0,1,ii,jj]
#                          femIntegrals.zeroJacobian(jacIndex,
#                                                    jacobian)
#                  for ii in range(self.l2g['nFreeDOF'][right_eN_global]):
#                      for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                          jacIndex = self.csrRowIndeces[right_eN_global,ii] + self.csrColumnOffsets_eb[ebN,1,0,ii,jj]
#                          femIntegrals.zeroJacobian(jacIndex,
#                                                    jacobian)
#                      for  jj  in range(self.l2g['nFreeDOF'][right_eN_global]):
#                          jacIndex = self.csrRowIndeces[right_eN_global,ii] + self.csrColumnOffsets_eb[ebN,1,1,ii,jj]
#                          femIntegrals.zeroJacobian(jacIndex,
#                                                    jacobian)
#              for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#                  ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#                  eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  for  ii  in range(self.l2g['nFreeDOF'][eN_global]):
#                      for  jj  in range(self.l2g['nFreeDOF'][eN_global]):
#                          jacIndex = self.csrRowIndeces[eN_global,ii] + self.csrColumnOffsets_eb[ebN,0,0,ii,jj]
#                          femIntegrals.zeroJacobian(jacIndex,
#                                                    jacobian)
        #beginAssembly(jacobian)
        #self.initializeJacobian(jacobian)
#          jacobianTemp = SparseMat(self.nFreeDOF_global,self.nFreeDOF_global,self.nFreeDOF_global)
#          for eN in range(self.mesh.nElements_global):
#              localDOF = self.l2g[eN]['freeLocal']
#              freeGlobalDOF = self.l2g[eN]['freeGlobal']
#              for i,I in zip(localDOF,freeGlobalDOF):
#                  for j,J in zip(localDOF,freeGlobalDOF):
#                      jacobian[I,J] += self.elementJacobian[eN,i,j]
#          for eN in range(self.mesh.nElements_global):
#              for ii in range(self.l2g['nFreeDOF'][eN]):
#                  I = self.l2g['freeGlobal'][eN,ii]
#                  for jj in range(self.l2g['nFreeDOF'][eN]):
#                      J = self.l2g['freeGlobal'][eN,jj]
#                      jacobian[I,J] =0.0
#          for eN in range(self.mesh.nElements_global):
#              for ii in range(self.l2g['nFreeDOF'][eN]):
#                  i=self.l2g['freeLocal'][eN,ii]
#                  I=self.l2g['freeGlobal'][eN,ii]
#                  for jj in range(self.l2g['nFreeDOF'][eN]):
#                      j=self.l2g['freeLocal'][eN,jj]
#                      J=self.l2g['freeGlobal'][eN,jj]
#                      jacobian[I,J] += self.elementJacobian[eN,i,j]
#         import inspect
#         print inspect.getmembers(jacobian)
        femIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.mesh.nElements_global,
                                                                 self.nDOF_element, 
								 self.l2g['nFreeDOF'],
                                                                 self.l2g['freeLocal'],
                                                                 self.csrRowIndeces,
                                                                 self.csrColumnOffsets,
                                                                 self.elementJacobian,
                                                                 jacobian)
#          for eN in range(self.mesh.nElements_global):
#              for ii in range(self.l2g['nFreeDOF'][eN]):
#                  for jj in range(self.l2g['nFreeDOF'][eN]):
#                      jacobian.val[self.csrRowIdeces[eN,ii]+self.csrColumnOffsets[eN,ii,jj]] = 0.0
#          print jacobian
#         print `jacobian`
#          if self.numericalFlux != None:
#              for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
#                  ebN = self.mesh.interiorElementBoundariesArray[ebNI]
#                  left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
#                  left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
#                  left_localDOF = self.l2g[left_eN_global]['freeLocal']
#                  left_freeGlobalDOF = self.l2g[left_eN_global]['freeGlobal']
#                  right_localDOF = self.l2g[right_eN_global]['freeLocal']
#                  right_freeGlobalDOF = self.l2g[right_eN_global]['freeGlobal']
#                  for i,left_I in zip(left_localDOF,left_freeGlobalDOF):
#                      for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                          for j,left_J in zip(left_localDOF,left_freeGlobalDOF):
#                              jacobian[left_I,left_J]   += self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,0,k,j]*self.ebq['w'][left_eN_global,left_ebN_element,k,i]
#                          for j,right_J in zip(right_localDOF,right_freeGlobalDOF):
#                              jacobian[left_I,right_J]  += self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,1,k,j]*self.ebq['w'][left_eN_global,left_ebN_element,k,i]
#                  for i,right_I in zip(right_localDOF,right_freeGlobalDOF):
#                      for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                          for j,left_J in zip(left_localDOF,left_freeGlobalDOF):
#                              jacobian[right_I,left_J]  -= self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,0,k,j]*self.ebq['w'][right_eN_global,right_ebN_element,k,i]
#                          for j,right_J in zip(right_localDOF,right_freeGlobalDOF):
#                              jacobian[right_I,right_J] -= self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,1,k,j]*self.ebq['w'][right_eN_global,right_ebN_element,k,i]
#              for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#                  ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#                  eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  localDOF = self.l2g[eN_global]['freeLocal']
#                  freeGlobalDOF = self.l2g[eN_global]['freeGlobal']
#                  for i,I in zip(localDOF,freeGlobalDOF):
#                      for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                          for j,J in zip(localDOF,freeGlobalDOF):
#                              jacobian[I,J] += self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,0,k,j]*self.ebq['w'][eN_global,ebN_element,k,i]
#         if self.numericalFlux != None:
#             femIntegrals.updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(self.mesh.nInteriorElementBoundaries_global,
#                                                                                          self.nDOF_element,
#                                                                                          self.nElementBoundaryQuadraturePoints_elementBoundary,
#                                                                                          self.mesh.nElementBoundaries_element,
#                                                                                          self.mesh.interiorElementBoundariesArray,
#                                                                                          self.mesh.elementBoundaryElementsArray,
#                                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                                          self.l2g['nFreeDOF'],
#                                                                                          self.l2g['freeLocal'],
#                                                                                          self.csrRowIndeces,
#                                                                                          self.csrColumnOffsets_eb,
#                                                                                          self.fluxJacobian,
# 											 self.ebq['w'],
#                                                                                          jacobian)
#             femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(self.mesh.nExteriorElementBoundaries_global,
#                                                                                          self.nDOF_element,
#                                                                                          self.nElementBoundaryQuadraturePoints_elementBoundary,
#                                                                                          self.mesh.nElementBoundaries_element,
#                                                                                          self.mesh.exteriorElementBoundariesArray,
#                                                                                          self.mesh.elementBoundaryElementsArray,
#                                                                                          self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                                          self.l2g['nFreeDOF'],
#                                                                                          self.l2g['freeLocal'],
#                                                                                          self.csrRowIndeces,
#                                                                                          self.csrColumnOffsets_eb,
#                                                                                          self.fluxJacobian,
#                                                                                          self.ebq['w'],
#                                                                                          jacobian)
#         elif self.fluxBoundaryConditions == 'outFlow':
#            femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(self.mesh.nExteriorElementBoundaries_global,
#                                                                                         self.nDOF_element,
#                                                                                         self.nElementBoundaryQuadraturePoints_elementBoundary,
#                                                                                         self.mesh.nElementBoundaries_element,
#                                                                                         self.mesh.exteriorElementBoundariesArray,
#                                                                                         self.mesh.elementBoundaryElementsArray,
#                                                                                         self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                                         self.l2g['nFreeDOF'],
#                                                                                         self.l2g['freeLocal'],
#                                                                                         self.csrRowIndeces,
#                                                                                         self.csrColumnOffsets_eb,
#                                                                                         self.fluxJacobian,
#                                                                                         self.ebq['w'],
#                                                                                         jacobian)
#              for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
#                  ebN = self.mesh.interiorElementBoundariesArray[ebNI]
#                  left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
#                  left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
#                  for ii in range(self.l2g['nFreeDOF'][left_eN_global]):
#                      i = self.l2g['freeLocal'][left_eN_global,ii]
#                      left_I = self.l2g['freeGlobal'][left_eN_global,ii]
#                      for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                          for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                              j = self.l2g['freeLocal'][left_eN_global,jj]
#                              left_J = self.l2g['freeGlobal'][left_eN_global,jj]
#                              jacIndex = self.csrRowIndeces[left_eN_global,ii] + self.csrColumnOffsets_eb[ebN,0,0,ii,jj]
#                              femIntegrals.updateAddJacobian(jacIndex,
#                                                            self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,0,k,j]*self.ebq['w'][left_eN_global,left_ebN_element,k,i],
#                                                            jacobian)
#                          for jj in range(self.l2g['nFreeDOF'][right_eN_global]):
#                              j = self.l2g['freeLocal'][right_eN_global,jj]
#                              right_J = self.l2g['freeGlobal'][right_eN_global,jj]
#                              jacIndex = self.csrRowIndeces[left_eN_global,ii] + self.csrColumnOffsets_eb[ebN,0,1,ii,jj]
#                              femIntegrals.updateAddJacobian(jacIndex,
#                                                            self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,1,k,j]*self.ebq['w'][left_eN_global,left_ebN_element,k,i],
#                                                            jacobian)
#                  for ii in range(self.l2g['nFreeDOF'][right_eN_global]):
#                      i = self.l2g['freeLocal'][right_eN_global,ii]
#                      right_I = self.l2g['freeGlobal'][right_eN_global,ii]
#                      for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                          for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                              j = self.l2g['freeLocal'][left_eN_global,jj]
#                              left_J = self.l2g['freeGlobal'][left_eN_global,jj]
#                              jacIndex = self.csrRowIndeces[right_eN_global,ii] + self.csrColumnOffsets_eb[ebN,1,0,ii,jj]
#                              femIntegrals.updateAddJacobian(jacIndex,
#                                                             -self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,0,k,j]*self.ebq['w'][right_eN_global,right_ebN_element,k,i],
#                                                             jacobian)
#                          for  jj  in range(self.l2g['nFreeDOF'][right_eN_global]):
#                              j = self.l2g['freeLocal'][right_eN_global,jj]
#                              right_J = self.l2g['freeGlobal'][right_eN_global,jj]
#                              jacIndex = self.csrRowIndeces[right_eN_global,ii] + self.csrColumnOffsets_eb[ebN,1,1,ii,jj]
#                              femIntegrals.updateAddJacobian(jacIndex,
#                                                             -self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,1,k,j]*self.ebq['w'][right_eN_global,right_ebN_element,k,i],
#                                                             jacobian)
#              for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#                  ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#                  eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  for  ii  in range(self.l2g['nFreeDOF'][eN_global]):
#                      i = self.l2g['freeLocal'][eN_global,ii]
#                      I = self.l2g['freeGlobal'][eN_global,ii]
#                      for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                          for  jj  in range(self.l2g['nFreeDOF'][eN_global]):
#                              j = self.l2g['freeLocal'][eN_global,jj]
#                              J = self.l2g['freeGlobal'][eN_global,jj]
#                              jacIndex = self.csrRowIndeces[eN_global,ii] + self.csrColumnOffsets_eb[ebN,0,0,ii,jj]
#                              femIntegrals.updateAddJacobian(jacIndex,
#                                                            self.ebq_global['dElementBoundaryFlux_dv*dx_f'][ebN,0,k,j]*self.ebq['w'][eN_global,ebN_element,k,i],
#                                                            jacobian)
        #endAssembly(jacobian)
        return jacobian
    def updateElementResidual(self):
        """Calculate all the element residuals"""
        self.elementResidual.flat[:]=0.0
        if self.coefficients.mass != None:
            self.updateMass(self.elementResidual)
        self.updateHamiltonian(self.elementResidual)
        if self.stabilization != None:
            self.updateStabilization(self.elementResidual)
        if self.shockCapturing != None:
            self.updateShockCapturing(self.elementResidual)
        if self.numericalFlux != None:
            if self.coefficients.advection != None:
                self.updateAdvectionInteriorElementBoundaryFlux(self.elementResidual)
            if self.coefficients.diffusion != None:
                self.updateDiffusionElementBoundaryFlux(self.elementResidual)
        if self.dirichletNodeSetList != None:
            for eN in range(self.mesh.nElements_global):
                for j in self.dirichletNodeSetList[eN]:
                    J = self.trialSpace.dofMap.l2g[eN,j]
                    self.u.dof[J] = self.dirichletValues[(eN,j)]
                    self.elementResidual[eN,j] = self.u.dof[J]-self.dirichletValues[(eN,j)]
    def updateElementJacobian(self):
        self.elementJacobian.flat[:]=0.0
        if self.coefficients.mass != None:
            if self.timeIntegration.massIsImplicit:
                self.updateMassJacobian(self.elementJacobian)
        if self.timeIntegration.hamiltonianIsImplicit:
            self.updateHamiltonianJacobian(self.elementJacobian)
        if self.stabilization != None:
            if self.timeIntegration.stabilizationIsImplicit:
                self.updateStabilizationJacobian(self.elementJacobian)
        if self.shockCapturing != None:
            if self.timeIntegration.shockCapturingIsImplicit:
                self.updateShockCapturingJacobian(self.elementJacobian)
        if self.timeIntegration.duStar_du != None:
            self.elementJacobian *= self.timeIntegration.duStar_du
        if self.dirichletNodeSetList != None:
            for eN in range(self.mesh.nElements_global):
                for j in self.dirichletNodeSetList[eN]:
                    self.elementJacobian[eN,j,:]=0.0
                    self.elementJacobian[eN,j,j]=1.0
    def updateCoefficients(self):
        self.updateElementCoefficients()
        if self.elementBoundaryIntegrals:
            self.updateElementBoundaryCoefficients()
    def updateElementCoefficients(self):
        """
        update the nonlinear coefficients at the quadrature points and nodes
        """
        #
        #get u and grad(u) at the quadrature points
        #
        self.calculateFiniteElementU()
        #self.u.getValues(self.quadraturePoints,self.q['u'])
        #self.u.getGradients(self.quadraturePoints,self.q['grad(u)'])
        #
        #get coefficients at the quadrature points
        #
        self.coefficients.evaluate(t  = self.T,
                                   x  = Numeric.reshape(self.q['x'].flat,(self.nQuadraturePoints_global,3)),
                                   u  = self.q['u'].flat,
                                   grad_u  = Numeric.reshape(self.q['grad(u)'].flat,(self.nQuadraturePoints_global,self.nSpace_global)),
                                   m  = self.q['m'].flat,
                                   dm = self.q['dm'].flat,
                                   h  = Numeric.reshape(self.q['h'].flat,(self.nQuadraturePoints_global,)),
                                   dh = Numeric.reshape(self.q['dh'].flat,(self.nQuadraturePoints_global,self.nSpace_global)),
                                   rh  = self.q['rh'].flat)
        #
        # m_t at the quadrature points
        #
        self.timeIntegration.updateMass(self.q['m'],self.q['mt'],self.q['dm'],self.q['dmt'])
        #
        # stabilization
        #
        if self.stabilization != None:
            femIntegrals.calculateStabilizationHJ(self.mesh.nElements_global,
                                                  self.nQuadraturePoints_element,
                                                  self.nSpace_global,
                                                  self.stabilization,
                                                  self.e['diameter'],
                                                  self.q['dh'],
                                                  self.q['dmt'],
                                                  self.q['cfl'],
                                                  self.q['tau'])
            self.calculateAdjoint()
            self.timeIntegration.updateStabilization(self.q['tau'])
            self.timeIntegration.updateAdjoint(self.q['Lstar*w'])
        else:
            self.q['cfl'].flat[:]=10
        #
        # L^* w at the quadrature points
        #
        self.timeIntegration.updateHamiltonian(self.q['h'],self.q['dh'],self.q['dhdu'],self.q['grad(u)'])
        self.timeIntegration.updateGradients(self.q['grad(u)*grad(w)'])
        femIntegrals.calculatePDEResidualHJ(self.mesh.nElements_global,
                                            self.nQuadraturePoints_element,
                                            self.nSpace_global,
                                            self.q['dh'],
                                            self.q['grad(u)'],
                                            self.q['rh'],
                                            self.q['mt'],
                                            self.q['pdeResidual'])
        femIntegrals.calculatePDEResidualJacobianHJ(self.mesh.nElements_global,
                                                    self.nQuadraturePoints_element,
                                                    self.nDOF_element,
                                                    self.nSpace_global,
                                                    self.q['dhdu'],
                                                    self.q['grad(v)'],
                                                    self.q['dmt'],
                                                    self.q['v'],
                                                    self.q['dpdeResidual'])
        if self.shockCapturing != None:
            femIntegrals.calculateShockCapturingHJ(self.mesh.nElements_global,
                                                   self.nQuadraturePoints_element,
                                                   self.nSpace_global,
                                                   self.shockCapturing,
                                                   self.shockCapturingDiffusion,
                                                   self.e['diameter'],
                                                   self.q['pdeResidual'],
                                                   self.q['mt'],
                                                   self.q['dh'],
                                                   self.q['numDiff'])
            self.timeIntegration.updateShockCapturing(self.q['numDiff'])
        if self.coefficients.mass != None:
            femIntegrals.calculateScalarScalarProduct(self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.q['mt'],
                                                      self.q['dx_m'],
                                                      self.q['mt*dx_m'])
            if self.timeIntegration.massIsImplicit:
                femIntegrals.calculateScalarScalarProduct(self.mesh.nElements_global,
                                                          self.nQuadraturePoints_element,
                                                          self.q['dmt'],
                                                          self.q['dx_m'],
                                                          self.q['dmt*dx_m'])
        femIntegrals.calculateScalarScalarProduct(self.mesh.nElements_global,
                                                  self.nQuadraturePoints_element,
                                                  self.q['h'],
                                                  self.q['dx_h'],
                                                  self.q['h*dx_h'])
        if self.timeIntegration.hamiltonianIsImplicit:
            femIntegrals.calculateVectorScalarProduct(self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.nSpace_global,
                                                      self.q['dh'],
                                                      self.q['dx_h'],
                                                      self.q['dh*dx_h'])
        if self.stabilization != None:
            femIntegrals.calculateScalarScalarProduct(self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.q['tau'],
                                                      self.q['dx_stab'],
                                                      self.q['tau*dx_stab'])
        if self.shockCapturing != None:
            femIntegrals.calculateScalarScalarProduct(self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.q['numDiff'],
                                                      self.q['dx_numDiff'],
                                                      self.q['numDiff*dx_numDiff'])
    def updateElementBoundaryCoefficients(self):
        """
        Update the nonlinear coefficients at the element boundary quadrature points
        """
        #
        #get u and grad(u) at the quadrature points
        #
        self.calculateFiniteElementU_elementBoundaries()
    def updateQuadrature(self):
        self.updateElementQuadrature()
        if self.elementBoundaryIntegrals:
            self.updateElementBoundaryQuadrature()
    
    def updateElementQuadrature(self):
        """
        Update the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.
        
        This function should be called only when the mesh changes.
        """
        #TODO: Extend to allow phi to be in a different space than v
        #
        #get physical locations of quadrature points and jacobian information there
        self.u.femSpace.elementMaps.getValues(self.quadraturePoints,
                                              self.q['x'])
        self.u.femSpace.elementMaps.getJacobianValues(self.quadraturePoints,
                                                      self.q['J'],
                                                      self.q['inverse(J)'],
                                                      self.q['det(J)'])
        #
        #get shape information at the quadrature points
        #
        self.testSpace.getBasisValues(self.quadraturePoints,
                                      self.q['w'])
        self.testSpace.getBasisGradientValues(self.quadraturePoints,
                                              self.q['inverse(J)'],
                                              self.q['grad(w)'])
        self.trialSpace.getBasisValues(self.quadraturePoints,
                                       self.q['v'])
        self.trialSpace.getBasisGradientValues(self.quadraturePoints,
                                               self.q['inverse(J)'],
                                               self.q['grad(v)'])
        #
        #get tensor products of shape information
        #
        if self.testCalculateShape:
            for eN in range(self.mesh.nElements_global):
                for k,xi in enumerate(self.quadraturePoints):
                    for i in range(self.nDOF_element):
                        for j in range(self.nDOF_element):
                            self.q['v*w'][eN,k,j,i] = self.q['v'][eN,k,j]*\
                                                         self.q['w'][eN,k,i]
                            for l in range(self.nSpace_global):
                                self.q['v*grad(w)'][eN,k,j,i,l] = self.q['v'][eN,k,j]*self.q['grad(w)'][eN,k,i,l]
                                #unused storage?
                                self.q['grad(v)*w'][eN,k,j,i,l] = self.q['grad(v)'][eN,k,j,l]*self.q['w'][eN,k,i]
                                for m in range(self.nSpace_global):
                                    self.q['grad(v)*grad(w)'][eN,k,j,i,l,m] = self.q['grad(v)'][eN,k,j,l]*self.q['grad(w)'][eN,k,i,m]
#              print self.q['v*w']
#              print self.q['v*grad(w)']
#         print "grad(v)*w old"
#         print self.q['grad(v)*w']
#         tmp = self.q['grad(v)*w'][:]
#              print self.q['grad(v)*grad(w)']
        femIntegrals.calculateShape_x_Shape(self.mesh.nElements_global,
                                            self.nQuadraturePoints_element,
                                            self.nDOF_element,
                                            self.q['v'],
                                            self.q['w'],
                                            self.q['v*w'])
        femIntegrals.calculateShape_x_GradShape(self.mesh.nElements_global,
                                                self.nQuadraturePoints_element,
                                                self.nDOF_element,
                                                self.nSpace_global,
                                                self.q['v'],
                                                self.q['grad(w)'],
                                                self.q['v*grad(w)'])
        femIntegrals.calculateGradShape_x_Shape(self.mesh.nElements_global,
                                                self.nQuadraturePoints_element,
                                                self.nDOF_element,
                                                self.nSpace_global,
                                                self.q['grad(v)'],
                                                self.q['w'],
                                                self.q['grad(v)*w'])
        femIntegrals.calculateGradShape_x_GradShape(self.mesh.nElements_global,
                                                    self.nQuadraturePoints_element,
                                                    self.nDOF_element,
                                                    self.nSpace_global,
                                                    self.q['grad(v)'],
                                                    self.q['grad(w)'],
                                                    self.q['grad(v)*grad(w)'])
        if self.testCalculateShape:
            print "from C code"
#              print self.q['v*w']
#              print self.q['v*grad(w)']
#              print self.q['grad(v)*w']
#              print self.q['grad(v)*grad(w)']
        #
        #get quadrature weights for each integral
        #
        self.q['abs(det(J))']=Numeric.absolute(self.q['det(J)'])
        for integral in self.integralKeys:
            femIntegrals.calculateIntegrationWeights(self.mesh.nElements_global,
                                                     self.nQuadraturePoints_element,
                                                     self.q['abs(det(J))'],
                                                     self.quadratureWeights[integral],
                                                     self.q['dx_'+integral])
    def updateElementBoundaryQuadrature(self):
        """
        Update the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on element boundaries.

        This function should be called only when the mesh changes.
        """
        #
        #get physical locations of element boundary quadrature points
        #
        self.u.femSpace.elementMaps.getValuesTrace(self.elementBoundaryQuadraturePoints,
                                                 self.ebq['x'])
        #
        #get metric tensor and unit normals
        #
        self.u.femSpace.elementMaps.getJacobianValuesTrace(self.elementBoundaryQuadraturePoints,
                                                         self.ebq['inverse(J)'],
                                                         self.ebq['g'],
                                                         self.ebq['sqrt(det(g))'],
                                                         self.ebq['n'])
        #TODO: vectorize these loops
        #
        #the boundary quadrature points will correspond to the ones on the "left side"
        #of each element boundary so we need to fix the "right side"
        #
        #first copy left information into ebq_global storage 
        for ebN in range(self.mesh.nElementBoundaries_global):
            left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                self.ebq_global['x'][ebN,k] = self.ebq['x'][left_eN_global,left_ebN_element,k]
                self.ebq_global['n'][ebN,k,:] = self.ebq['n'][left_eN_global,left_ebN_element,k]
        #now copy left information into right physical points on interior element boundaries
        for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
            ebN = self.mesh.interiorElementBoundariesArray[ebNI]
            right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
            right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                self.ebq['x'][right_eN_global,right_ebN_element,k,:] = self.ebq_global['x'][ebN,k]
        #now map the physical points back to the reference element
        self.u.femSpace.elementMaps.getInverseValuesTrace(self.ebq['x'],self.ebq['hat(x)'])
        #
        #since the points on the reference boundary may be reordered on many right element boundaries, we
        #have to use an array of reference boundary points on all element boundaries
        #first copy the left reference element boundary quadrature points from the reference element boundary
        for ebN in range(self.mesh.nElementBoundaries_global):
            left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                self.ebq['bar(x)'][left_eN_global,left_ebN_element,k,:] = self.elementBoundaryQuadraturePoints[k]
        #now get the right element boundary quadrature points and normals from the inverse map
        #note that the inverses are only defined on the boundaries of the reference element
        boundaryMapInverseList = self.u.femSpace.referenceFiniteElement.referenceElement.boundaryMapInverseList
        for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
            ebN = self.mesh.interiorElementBoundariesArray[ebNI]
            left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
            left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                self.ebq['bar(x)'][right_eN_global,right_ebN_element,k,:] = boundaryMapInverseList[right_ebN_element](
                    self.ebq['hat(x)'][right_eN_global,right_ebN_element,k])
                self.ebq['n'][right_eN_global,right_ebN_element,k] = - self.ebq['n'][left_eN_global,left_ebN_element,k]
        #
        #get the shape information at the reference element boundary quadrature points 
        #   
        self.testSpace.getBasisValuesTraceAtArray(self.ebq['bar(x)'],
                                                  self.ebq['w'])
        self.trialSpace.getBasisValuesTraceAtArray(self.ebq['bar(x)'],
                                                   self.ebq['v'])
        self.trialSpace.getBasisGradientValuesTraceAtArray(self.ebq['bar(x)'],
                                                           self.ebq['inverse(J)'],
                                                           self.ebq['grad(v)'])
    def calculateAdjoint(self):
        """
        Calculate the action of the adjoint of the linearized spatial
        operator on the test functions.
        """
        femIntegrals.calculateAdjointHJ(self.mesh.nElements_global,
                                        self.nQuadraturePoints_element,
                                        self.nDOF_element,
                                        self.nSpace_global,
                                        self.q['w'],
                                        self.q['grad(w)'],
                                        self.q['dh'],
                                        self.q['rh'],
                                        self.q['Lstar*w'])
    def calculateFiniteElementU(self):
        """
        Calculate u, grad(u), and the (outer) product grad(u) x grad(w)
        """
        self.u.getValues(self.q['v'],
                         self.q['u'])
        self.u.getGradientValues(self.q['grad(v)'],
                                 self.q['grad(u)'])
        if  self.shockCapturing != None:
            self.u.getGradientTensorValues(self.q['grad(v)*grad(w)'],
                                           self.q['grad(u)*grad(w)'])
    def calculateFiniteElementU_elementBoundaries(self):
        """
        Calculate u and grad(u)
        """
        self.u.getValuesTrace(self.ebq['v'],self.ebq['u'])
        self.u.getGradientValuesTrace(self.ebq['grad(v)'],self.ebq['grad(u)'])
    def updateZeroLevelSetDirichletConditions(self):
        self.dirichletNodeSetList = []
        self.dirichletGlobalNodeSet=set()
        self.dirichletValues={}
        for eN in range(self.mesh.nElements_global):
            self.dirichletNodeSetList.append(set())
            signU = 0
            j0=0
            while ((signU == 0) and
                   (j0 < self.nDOF_element)):
                J0 = self.trialSpace.dofMap.l2g[eN,j0]
                if self.u.dof[J0] < 0.0:
                    signU = -1
                elif  self.u.dof[J0] > 0.0:
                    signU = 1
                else:
                    self.dirichletNodeSetList[eN].add(j0)
                    self.dirichletValues[(eN,j0)]=float(self.u.dof[J0])
                    self.dirichletGlobalNodeSet.add(J0)
                j0 += 1
            for j in range(j0,self.nDOF_element):
                J = self.trialSpace.dofMap.l2g[eN,j]
                if (((self.u.dof[J] < 0.0) and
                     (signU == 1)) or
                    ((self.u.dof[J] > 0.0) and
                     (signU == -1))):
                    for jj in range(self.nDOF_element):
                        JJ = self.trialSpace.dofMap.l2g[eN,jj]
                        self.dirichletNodeSetList[eN].add(jj)
                        self.dirichletValues[(eN,jj)]=float(self.u.dof[JJ])
                        self.dirichletGlobalNodeSet.add(JJ)
                    break
                elif (self.u.dof[J] == 0.0):
                    self.dirichletNodeSetList[eN].add(j)
                    self.dirichletValues[(eN,j)]=float(self.u.dof[J])
                    self.dirichletGlobalNodeSet.add(J)
        for eN in range(self.mesh.nElements_global):
            for j in range(self.nDOF_element):
                J = self.trialSpace.dofMap.l2g[eN,j]
                if J in self.dirichletGlobalNodeSet:
                    self.dirichletNodeSetList[eN].add(j)
                    self.dirichletValues[(eN,j)]=float(self.u.dof[J])
    def updateMass(self,residual):
        """
        Calculate the element mass accumulation integral and update
        the element residual
        """
        if self.testUpdateIntegrals:
            self.restemp[:]=residual
        femIntegrals.updateMass(self.mesh.nElements_global,
                                self.nQuadraturePoints_element,
                                self.nDOF_element,
                                self.q['mt*dx_m'],
                                self.q['w'],
                                residual)
        if self.testUpdateIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        self.restemp[eN,i] += (self.q['mt*dx_m'][eN,k]*
                                           self.q['w'][eN,k,i])
            if min(self.restemp.flat == residual.flat) == 0:
                print "mass error"
                print self.restemp - residual
                print "mass error"
    def updateMassJacobian(self,jacobian):
        """
        Calculate the jacobian of the element mass accumulation
        integral and update the element jacobian
        """
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp.flat[:]=jacobian.flat
        femIntegrals.updateMassJacobian(self.mesh.nElements_global,
                                        self.nQuadraturePoints_element,
                                        self.nDOF_element,
                                        self.q['dmt*dx_m'],
                                        self.q['v*w'],
                                        jacobian)
        if self.testUpdateJacobianIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for j in range(self.nDOF_element):
                            self.jacobiantemp[eN,i,j] += (self.q['dmt*dx_m'][eN,k]*
                                                          self.q['v*w'][eN,k,j,i])
            if min(self.jacobiantemp.flat == jacobian.flat) == 0:
                print "mass jacobian error"
                print self.jacobiantemp - jacobian
                print "mass jacobian error"
    def updateHamiltonian(self,residual):
        """
        Calculate the element hamiltonian integral and update element
        residual
        """
        femIntegrals.updateHamiltonian(self.mesh.nElements_global,
                                       self.nQuadraturePoints_element,
                                       self.nDOF_element,
                                       self.q['h*dx_h'],
                                       self.q['w'],
                                       residual)
    def updateHamiltonianJacobian(self,jacobian):
        """
        Calculate the jacobian of the element advection integral and
        update element jacobian
        """
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp[:]=jacobian
        femIntegrals.updateHamiltonianJacobian(self.mesh.nElements_global,
                                               self.nQuadraturePoints_element,
                                               self.nDOF_element,
                                               self.nSpace_global,
                                               self.q['dh*dx_h'],
                                               self.q['grad(v)*w'],
                                               jacobian)
        if self.testUpdateJacobianIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for j in range(self.nDOF_element):
                            self.jacobiantemp[eN][i][j] += (
                                Numeric.dot(self.q['dh*dx_h'][eN][k],
                                            self.q['grad(v)*w'][eN][k][j][i]))
            if min(self.jacobiantemp.flat == jacobian.flat) == 0:
                print "hamiltonian jacobian error"
                print self.jacobiantemp - jacobian
                print "hamiltonian jacobian error"
    def updateStabilization(self,residual):
        """
        Calculate the element stabilization integral and update the element
        residual
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateStabilization(self.mesh.nElements_global,
                                         self.nQuadraturePoints_element,
                                         self.nDOF_element,
                                         self.q['tau*dx_stab'],
                                         self.q['pdeResidual'],
                                         self.q['Lstar*w'],
                                         residual)
        if self.testUpdateIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        self.restemp[eN,i] -= (self.q['tau*dx_stab'][eN,k]*
                                            self.q['pdeResidual'][eN,k]*
                                            self.q['Lstar*w'][eN,k,i])
            if min(self.restemp.flat == residual.flat) == 0:
                print "stabilization error"
                print self.restemp - residual
                print "stabilization error"
    def updateStabilizationJacobian(self,jacobian):
        """
        Calculate the jacobian of the element stabilization integral and
        update element jacobian
        """
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp.flat[:]=jacobian.flat
        femIntegrals.updateStabilizationJacobian(self.mesh.nElements_global,
                                                 self.nQuadraturePoints_element,
                                                 self.nDOF_element,
                                                 self.q['tau*dx_stab'],
                                                 self.q['dpdeResidual'],
                                                 self.q['Lstar*w'],
                                                 jacobian)
        if self.testUpdateJacobianIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for j in range(self.nDOF_element):
                            self.jacobiantemp[eN,i,j] -= (self.q['tau*dx_stab'][eN,k]*
                                                 self.q['dpdeResidual'][eN,k,j]*
                                                 self.q['Lstar*w'][eN,k,i])
            if min(self.jacobiantemp.flat == jacobian.flat) == 0:
                print "stabilization jacobian error"
                print self.jacobiantemp - jacobian
                print "stabilizatioon jacobian error"
    def updateShockCapturing(self,residual):
        """
        Calculate the element shock capturing integral and the update
        element residual
        """
        if self.testUpdateIntegrals:
            self.restemp[:]=residual
        femIntegrals.updateShockCapturing(self.mesh.nElements_global,
                                          self.nQuadraturePoints_element,
                                          self.nDOF_element,
                                          self.nSpace_global,
                                          self.q['numDiff*dx_numDiff'],
                                          self.q['grad(u)*grad(w)'],
                                          residual)
        if self.testUpdateIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for l in range(self.nSpace_global):
                            self.restemp[eN,i] += (self.q['numDiff*dx_numDiff'][eN,k]*
                                                self.q['grad(u)*grad(w)'][eN,k,i,l,l])
            if min(self.restemp.flat == residual.flat) == 0:
                print "shock capturing error"
                print self.restemp - residual
                print "shock capturing error"
    def updateShockCapturingJacobian(self,jacobian):
        """
        Calculate the jacobian of the element shock capturing integral
        and update the element jacobian
        """
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp.flat[:]=jacobian.flat
        femIntegrals.updateShockCapturingJacobian(self.mesh.nElements_global,
                                                  self.nQuadraturePoints_element,
                                                  self.nDOF_element,
                                                  self.nSpace_global,
                                                  self.q['numDiff*dx_numDiff'],
                                                  self.q['grad(v)*grad(w)'],
                                                  jacobian)
        if self.testUpdateJacobianIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for j in range(self.nDOF_element_phi):
                            J = self.trialSpace.dofMap.l2g[eN,j]
                            for l in range(self.nSpace_global):
                                self.jacobiantemp[eN,i,j] += (self.q['numDiff*dx_numDiff'][eN,k]*
                                                     self.q['grad(v)*grad(w)'][eN,k,j,i,l,l])
            if min(self.jacobiantemp.flat == jacobian.flat) == 0:
                print "shock capturing jacobian error"
                print self.jacobiantemp - jacobian
                print "shock capturing jacobian error"
    def setFreeDOF(self,free_u):
        """
        Set the free(non-Dirichlet) DOF from the global DOF
        """
        for free_dofN,dofN in enumerate(self.dirichletConditions.freeDOFSet):
            free_u[free_dofN] = self.u.dof[dofN]
    def updateLocal2Global(self):
        """
        Build a mapping between local DOF numbers and global free DOF
        numbers, that is, EXCLUDING Dirichlet DOF.  We have to use a
        list of local dof because of the exclusion of Dirichlet nodes
        (otherwise we could just loop over range(self.nDOF_element).
        """
        self.l2g={'nFreeDOF':Numeric.zeros((self.mesh.nElements_global,),Numeric.Int),
                  'freeLocal':Numeric.zeros((self.mesh.nElements_global,self.nDOF_element),Numeric.Int),
                  'freeGlobal':Numeric.zeros((self.mesh.nElements_global,self.nDOF_element),Numeric.Int)}
#          for eN in range(self.mesh.nElements_global):
#              self.l2g[eN]={'freeLocal':[],'freeGlobal':[]}
#              for i in range(self.nDOF_element):
#                  I = self.testSpace.dofMap.l2g[eN,i]
#                  if self.dirichletConditions.global2freeGlobal.has_key(I):
#                      self.l2g[eN]['freeLocal'].append(i)
#                      self.l2g[eN]['freeGlobal'].append(self.dirichletConditions.global2freeGlobal[I])
        for eN in range(self.mesh.nElements_global):
            nFreeDOF=0
            for i in range(self.nDOF_element):
                I = self.testSpace.dofMap.l2g[eN,i]
                if self.dirichletConditions.global2freeGlobal.has_key(I):
                    self.l2g['freeLocal'][eN,nFreeDOF]=i
                    self.l2g['freeGlobal'][eN,nFreeDOF]=self.dirichletConditions.global2freeGlobal[I]
                    nFreeDOF+=1
            self.l2g['nFreeDOF'][eN]=nFreeDOF
    def setUnknowns(self,free_u):
        #Load the unknowns
        for free_dofN,dofN in enumerate(self.dirichletConditions.freeDOFSet):
            self.u.dof[dofN]=free_u[free_dofN]   
    def initializeJacobian(self):
#          for eN in range(self.mesh.nElements_global):
#              localDOF = self.l2g[eN]['freeLocal']
#              freeGlobalDOF = self.l2g[eN]['freeGlobal']
#              for i,I in zip(localDOF,freeGlobalDOF):
#                  for j,J in zip(localDOF,freeGlobalDOF):
#                      jacobian[I,J] =0.0
        #build sparse  matrix and extract the stuff we  need  to directly load/update csr matrices
        #from element matrices
        columnIndecesDict={}
        for eN in range(self.mesh.nElements_global):
            for ii in range(self.l2g['nFreeDOF'][eN]):
                I = self.l2g['freeGlobal'][eN,ii]
                if not columnIndecesDict.has_key(I):
                    columnIndecesDict[I]=set()
                for jj in range(self.l2g['nFreeDOF'][eN]):
                    J = self.l2g['freeGlobal'][eN,jj]
                    columnIndecesDict[I].add(J)
                    #jacobian[I,J] =1.0
#         if self.numericalFlux != None:
#             for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
#                 ebN = self.mesh.interiorElementBoundariesArray[ebNI]
#                 left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                 right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
#                 left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                 right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
#                 for ii in range(self.l2g['nFreeDOF'][left_eN_global]):
#                     i = self.l2g['freeLocal'][left_eN_global,ii]
#                     left_I = self.l2g['freeGlobal'][left_eN_global,ii]
#                     if not columnIndecesDict.has_key(left_I):
#                         columnIndecesDict[left_I]=set()
#                     for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                         j = self.l2g['freeLocal'][left_eN_global,jj]
#                         left_J = self.l2g['freeGlobal'][left_eN_global,jj]
#                         columnIndecesDict[left_I].add(left_J)
#                         jacobian[left_I,left_J]   =1.0
#                     for jj in range(self.l2g['nFreeDOF'][right_eN_global]):
#                         j = self.l2g['freeLocal'][right_eN_global,jj]
#                         right_J = self.l2g['freeGlobal'][right_eN_global,jj]
#                         columnIndecesDict[left_I].add(right_J)
#                         jacobian[left_I,right_J]  =1.0
#                 for ii in range(self.l2g['nFreeDOF'][right_eN_global]):
#                     i = self.l2g['freeLocal'][right_eN_global,ii]
#                     right_I = self.l2g['freeGlobal'][right_eN_global,ii]
#                     if not columnIndecesDict.has_key(right_I):
#                         columnIndecesDict[right_I]=set()
#                     for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                         j = self.l2g['freeLocal'][left_eN_global,jj]
#                         left_J = self.l2g['freeGlobal'][left_eN_global,jj]
#                         columnIndecesDict[right_I].add(left_J)
#                         jacobian[right_I,left_J]  =1.0
#                     for  jj  in range(self.l2g['nFreeDOF'][right_eN_global]):
#                         j = self.l2g['freeLocal'][right_eN_global,jj]
#                         right_J = self.l2g['freeGlobal'][right_eN_global,jj]
#                         columnIndecesDict[right_I].add(right_J)
#                         jacobian[right_I,right_J] =1.0
#         if self.numericalFlux != None  or  (self.fluxBoundaryConditions == 'outFlow'):
#             for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#                 ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#                 eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                 ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                 for  ii  in range(self.l2g['nFreeDOF'][eN_global]):
#                     i = self.l2g['freeLocal'][eN_global,ii]
#                     I = self.l2g['freeGlobal'][eN_global,ii]
#                     if not columnIndecesDict.has_key(I):
#                         columnIndecesDict[I]=set()
#                     for  jj  in range(self.l2g['nFreeDOF'][eN_global]):
#                         j = self.l2g['freeLocal'][eN_global,jj]
#                         J = self.l2g['freeGlobal'][eN_global,jj]
#                         columnIndecesDict[I].add(J)
#                         jacobian[I,J] =1.0
        nnz = 0
        rowptr = Numeric.zeros(self.nFreeDOF_global+1,Numeric.Int)
        lastIndex=0
        columnOffsetDict={}
        for I in range(self.nFreeDOF_global):
            columnIndeces=list(columnIndecesDict[I])
            columnIndeces.sort()
            rowptr[I]=lastIndex
            lastIndex += len(columnIndeces)
        rowptr[self.nFreeDOF_global]=lastIndex
        nnz = lastIndex
        colind = Numeric.zeros((nnz,),Numeric.Int)
        for I in range(self.nFreeDOF_global):
            columnIndeces=list(columnIndecesDict[I])
            for columnOffset,J in enumerate(columnIndeces):
                columnOffsetDict[(I,J)] = columnOffset
                colind[rowptr[I]+columnOffset]=J
        if self.matType == 'csr':
            self.nzval = Numeric.zeros((nnz,),Numeric.Float)
            self.jacobian = SparseMat(self.nFreeDOF_global,self.nFreeDOF_global,nnz,self.nzval,colind,rowptr)
        else:
            self.jacobian = Mat(self.nFreeDOF_global,self.nFreeDOF_global)
        self.csrRowIndeces = Numeric.zeros((self.mesh.nElements_global,self.nDOF_element),Numeric.Int)
        self.csrColumnOffsets = Numeric.zeros((self.mesh.nElements_global,self.nDOF_element,self.nDOF_element),Numeric.Int)
        self.csrColumnOffsets_eb = Numeric.zeros((self.mesh.nElementBoundaries_global,2,2,self.nDOF_element,self.nDOF_element),Numeric.Int)
        for eN in range(self.mesh.nElements_global):
            for ii in range(self.l2g['nFreeDOF'][eN]):
                I = self.l2g['freeGlobal'][eN,ii]
                self.csrRowIndeces[eN,ii]=rowptr[I]
                for jj in range(self.l2g['nFreeDOF'][eN]):
                    J = self.l2g['freeGlobal'][eN,jj]
                    self.csrColumnOffsets[eN,ii,jj] = columnOffsetDict[(I,J)]
#         if self.numericalFlux != None:
#             for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
#                 ebN = self.mesh.interiorElementBoundariesArray[ebNI]
#                 left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                 right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
#                 left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                 right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
#                 for ii in range(self.l2g['nFreeDOF'][left_eN_global]):
#                     i = self.l2g['freeLocal'][left_eN_global,ii]
#                     left_I = self.l2g['freeGlobal'][left_eN_global,ii]
#                     for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                         j = self.l2g['freeLocal'][left_eN_global,jj]
#                         left_J = self.l2g['freeGlobal'][left_eN_global,jj]
#                         self.csrColumnOffsets_eb[ebN,0,0,ii,jj] = columnOffsetDict[(left_I,left_J)]
#                     for jj in range(self.l2g['nFreeDOF'][right_eN_global]):
#                         j = self.l2g['freeLocal'][right_eN_global,jj]
#                         right_J = self.l2g['freeGlobal'][right_eN_global,jj]
#                         self.csrColumnOffsets_eb[ebN,0,1,ii,jj] = columnOffsetDict[(left_I,right_J)]
#                 for ii in range(self.l2g['nFreeDOF'][right_eN_global]):
#                     i = self.l2g['freeLocal'][right_eN_global,ii]
#                     right_I = self.l2g['freeGlobal'][right_eN_global,ii]
#                     for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
#                         j = self.l2g['freeLocal'][left_eN_global,jj]
#                         left_J = self.l2g['freeGlobal'][left_eN_global,jj]
#                         self.csrColumnOffsets_eb[ebN,1,0,ii,jj] = columnOffsetDict[(right_I,left_J)]
#                     for  jj  in range(self.l2g['nFreeDOF'][right_eN_global]):
#                         j = self.l2g['freeLocal'][right_eN_global,jj]
#                         right_J = self.l2g['freeGlobal'][right_eN_global,jj]
#                         self.csrColumnOffsets_eb[ebN,1,1,ii,jj] = columnOffsetDict[(right_I,right_J)]
#         if self.numericalFlux != None  or (self.fluxBoundaryConditions == 'outFlow'):
#             for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#                 ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#                 eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                 ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                 for  ii  in range(self.l2g['nFreeDOF'][eN_global]):
#                     i = self.l2g['freeLocal'][eN_global,ii]
#                     I = self.l2g['freeGlobal'][eN_global,ii]
#                     for  jj  in range(self.l2g['nFreeDOF'][eN_global]):
#                         j = self.l2g['freeLocal'][eN_global,jj]
#                         J = self.l2g['freeGlobal'][eN_global,jj]
#                         self.csrColumnOffsets_eb[ebN,0,0,ii,jj] = columnOffsetDict[(I,J)]
        self.nNonzerosInJacobian = len(columnOffsetDict)
        #print self.csrRowIndeces
        #print self.csrColumnOffsets
        #done building csr data structures
        return self.jacobian

class MultilevelHamiltonJacobi:
    """Nonlinear ADR on a multilevel mesh"""
    
    def __init__(self,nd,mlMesh,
                 TrialSpaceType,
                 TestSpaceType,
                 matType,
                 dirichletConditionSetter,
                 coefficients,
                 quadrature,
                 elementBoundaryQuadrature,
                 stabilization='2',
                 shockCapturing='1',
                 shockCapturingDiffusion=0.1,
                 TimeIntegrationClass = BackwardEuler,
                 calculateElementBoundaryValues = False,
		 tOrder=1):

        """read in the multilevel mesh, mesh independent boundary
        conditions, and types for test and trial spaces and the
        jacobian. Pass through the rest to the models on each mesh"""
        self.jacobianList=[]
        self.uList=[]
        self.duList=[]
        self.rList=[]
        self.trialSpaceList=[]
        self.bcList=[]
        self.modelList=[]
        print "Building HamiltonJacobi for each mesh"
        for mesh in mlMesh.meshList:
            print "Generating Trial Space"
            trialSpace = TrialSpaceType(mesh,nd)
            self.trialSpaceList.append(trialSpace)
            print "Generating Test Space"
            testSpace = TestSpaceType(mesh,nd)
            print "Allocating u"
            u = FiniteElementFunction(trialSpace)
            print "Setting Boundary Conditions"
            dirichletConditions=DOFBoundaryConditions(
                trialSpace,dirichletConditionSetter)
            self.bcList.append(dirichletConditions)
            print "Initializing OneLevelHamiltonJacobi"
            hamiltonJacobi=OneLevelHamiltonJacobi(u,
                                                  testSpace,
                                                  matType,
                                                  dirichletConditions,
                                                  copy.deepcopy(coefficients),
                                                  quadrature,
                                                  elementBoundaryQuadrature,
                                                  stabilization,
                                                  shockCapturing,
                                                  shockCapturingDiffusion,
                                                  TimeIntegrationClass,
                                                  calculateElementBoundaryValues,
						  tOrder)
            self.modelList.append(hamiltonJacobi)
            print "Allocating Vectors"
            self.u = Numeric.zeros((hamiltonJacobi.dim,),Numeric.Float)
            self.du = Numeric.zeros((hamiltonJacobi.dim,),Numeric.Float)
            self.r = Numeric.zeros((hamiltonJacobi.dim,),Numeric.Float)
            print "Allocating Jacobian"
            self.jacobian = hamiltonJacobi.initializeJacobian()
            self.jacobianList.append(self.jacobian)
            self.uList.append(self.u)
            self.duList.append(self.du)
            self.rList.append(self.r)
        print "Building Mesh Transfers"
        self. meshTransfers = MultilevelProjectionOperators(mlMesh,
							    TrialSpaceType,
							    nd,
							    self.bcList,
							    self.trialSpaceList)
    def setInitialConditions(self,getInitialConditions,T=0.0):
        for m,u in zip(self.modelList,self.uList):
            m.setInitialConditions(getInitialConditions,T)
            m.setFreeDOF(u)
    def initializeTimeIntegration(self):
        for m in self.modelList:
            m.initializeTimeIntegration()
    def chooseDT(self,DTSET=None):
        if DTSET == None:
            self.modelList[-1].timeIntegration.chooseDT()
        else:
            self.modelList[-1].timeIntegration.DT = DTSET
        self.DT = self.modelList[-1].timeIntegration.DT
        for m in self.modelList:
            m.timeIntegration.DT = self.DT
            m.T += self.DT
    def updateTimeHistory(self):
        for m,u,r in zip(self.modelList,self.uList,self.rList):
            m.timeIntegration.updateTimeHistory()

## @}

if __name__ == '__main__':
    import profile
    import pstats
    import HamiltonJacobiTests
    profile.run('HamiltonJacobiTests.runTests()','HamiltonJacobiProf')
    p = pstats.Stats('HamiltonJacobiProf')
    p.sort_stats('cumulative').print_stats(20)
    p.sort_stats('time').print_stats(20)
