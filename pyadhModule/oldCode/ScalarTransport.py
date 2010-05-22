import copy
from math import *
from EGeometry import *
from LinearAlgebra import *
from LinearSolvers import *
from MeshTools import *
from FemTools import *
from QuadTools import *
from TimeIntegrationTools import *
from NonlinearSolvers import *
import femIntegrals
## \defgroup ScalarTransport ScalarTransport
#
# @{

"""
This module contains methods for solving

m_t + \deld (f - a \grad \phi) + r = 0

The solution is u, and the nonlinear coefficients
are m,f,a,\phi, and r.
"""


class ScalarTransportCoefficients:
    """
    This the base class for coefficients of the scalar transport
    equation. The coefficients are evaluated at time t over a 1D array
    of points, x, and solution values, u. To create a numerical model for a scalar
    nonlinear transport equation it may only be necessary to derive a
    class from this base class and pass an instance to the appropriate
    scalar transport class.
    
    members:
    evaluate(...)
    
    string flags:
    mass      -- None,'linear', or 'nonlinear'
    advection -- None,'linear', or 'nonlinear'
    diffusion -- None,'constant', 'linear' or 'nonlinear'
    potential -- None,'constant', 'linear' or 'nonlinear'
    reaction  -- None,'linear', or 'nonlinear'
    
    These flags are used to allow optimizations in
    special cases.
    """
    def __init__(self,
                 mass='nonlinear',
                 advection='nonlinear',
                 diffusion='nonlinear',
                 potential='nonlinear',
                 reaction='nonlinear'):
        self.mass=mass
        self.advection=advection
        self.diffusion=diffusion
        self.potential=potential
        self.reaction=reaction
    def evaluate(self,
                 t,x,u,
                 m,dm,
                 f,df,
                 a,da,
                 phi,dphi,
                 r,dr):
        pass

class LinearADR_ConstantCoefficients(ScalarTransportCoefficients):
    """
    This class implements constant coefficients:

    (Mu)_t + \deld (Bu - A \grad u) + C u = 0 
    """
    from transportCoefficients import linearADR_ConstantCoefficientsEvaluate
    def __init__(self,M,A,B,C):
        ScalarTransportCoefficients.__init__(self,
                                             mass='linear',
                                             advection='linear',
                                             diffusion='constant',
                                             potential='linear',
                                             reaction='linear')
        self.M = M
        self.A = A
        self.B = B
        self.C = C
        self.useC=True#False
    def evaluate(self,
                 t,
                 x,
                 u,
                 m,dm,
                 f,df,
                 a,da,
                 phi,dphi,
                 r,dr):
        if self.useC:
            self.linearADR_ConstantCoefficientsEvaluate(x.shape[0],
                                                        f.shape[1],
                                                        self.M,
                                                        self.A,
                                                        self.B,
                                                        self.C,
                                                        t,
                                                        x,
                                                        u,
                                                        m,dm,
                                                        f,df,
                                                        a,da,
                                                        phi,dphi,
                                                        r,dr)
        else:
            m[:] = u
            m *= self.M
            dm[:] = self.M

            f[:]=u[:,Numeric.newaxis]
            f*=self.B
            df[:] = self.B

            a[:,:] = self.A
            da[:,:] = 0.0

            phi[:] = u
            dphi[:] = 1.0

            r[:] = u
            r*=self.C
            dr[:] = self.C

class OneLevelScalarTransport(NonlinearEquation):
    
    """
    A class for finite element discretizations of 
    advective-diffusive-reactive transport on a single spatial mesh.
    
    Objects of this type take the initial-boundary value
    problem for

    m_t + \deld (f - a \grad \phi) + r

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
 	
    """keys for each integral in the residual (element interior)"""
    integralKeys = ['m','f','a','r']
    
    """keys for each integral in the residual (element boundaries)"""
    elementBoundaryIntegralKeys = ['f','a']
    
    """keys for each coefficient of the pde (element interiors) """
    coefficientKeys = ['m','dm','f','df','a','da','phi','dphi','r','dr']
    
    """keys for each coefficient of the pde (element boundaries)"""
    elementBoundaryCoefficientKeys = ['f','df','a','da','phi','dphi']
    
    def __init__(self,
                 u,
                 phi,
                 testSpace,
                 matType,
                 dofBoundaryConditions,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditions='noFlow',#'outFlow','mixedFlow','setFlow'
                 stabilization='2',
                 shockCapturing=None,
                 shockCapturingDiffusion=0.1,
                 conservativeFlux=None,
                 numericalFlux=None,
                 TimeIntegrationClass=BackwardEuler,
                 tOrder = 1):

        """
        Allocate storage and initialize some variables.
        
        u   -- a FiniteElementFunction object
        
        phi -- a FiniteElmentFuncion object
        
        testSpace -- a FiniteElementSpace object, dim(testSpace) =
        dim(u.femSpace)
        
        dofBoundaryConditions -- a DOFBoundaryConditions object for
        the Dirichlet conditions
                
        coefficients -- a ScalarTransportCoefficients object
        
        elementQuadrature -- a dictionary of quadrature rules for each
        element integral
        
        elementBoundaryQuadrature -- a dictionary of quadrature rules
        for each element boundary integral
        
        stabilization -- an optional string specifying the stabilization method
        
        shockCapturing -- an optional string specifying the shock capturing method
        
        numericalFlux -- an optional string specifying the numerical flux calculation
        
        The constructor sets the input arguments, calculates some
        dimensions, and allocates some storage. The meanings of variable
        suffixes are
        
        _global          -- per physical domain
        _element         -- per element
        _elementBoundary -- per element boundary
        
        The prefix n means 'number of'.
        
        Storage is divided into quantities required at different sets
        of points or geometric entities. Each type of storage has a
        dictionary for all the quantities of that type. The names
        and dimensions of the storage dictionaries are
        
        e          -- at element
        eb_global  -- at element boundary, unique to element boundary
        q          -- at element quadrature        
        ebq        -- at element boundary quadrature, unique to elements
        ebq_global -- at element boundary quadrature, unique to element boundary
        phi_ip     -- at interpolation points for obtaing phi
        """
	
	#keys for each integral in the residual (element interior)"""
	self.integralKeys = ['m','f','a','r']
	
	#keys for each integral in the residual (element boundaries)"""
	self.elementBoundaryIntegralKeys = ['f','a']
	
	#keys for each coefficient of the pde (element interiors) """
	self.coefficientKeys = ['m','dm','f','df','a','da','phi','dphi','r','dr']
	
	#keys for each coefficient of the pde (element boundaries)"""
	self.elementBoundaryCoefficientKeys = ['f','df','a','da','phi','dphi']

        #
        #set the objects describing the method and boundary conditions
        #
        self.u = u
        self.phi = phi
        #Since phi is the interpolatory projection of phi(u) onto a finite
        #element space, we will need the derivative of the degrees of freedom
        #of phi with respect to the degrees of freedom of u. dPHI/dU.
        #For now, I'll assume that the finite element spaces
        #are the same for phi and u so that the dPHI/dU is simply
        #dphi_du(U), but I'm trying to work away from this assumption
        #so in some parts of the code things are more general. 
        #For now we abuse the FiniteElementFunction class for dPHI/dU
        self.dphi = FiniteElementFunction(phi.femSpace)
        self.trialSpace = self.u.femSpace
        self.matType = matType
        self.mesh = self.u.femSpace.elementMaps.mesh
        self.testSpace = testSpace
        self.dirichletConditions = dofBoundaryConditions
        self.dirichletNodeSetList=None
        self.coefficients = coefficients
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.shockCapturingDiffusion = shockCapturingDiffusion
        self.conservativeFlux = conservativeFlux
        self.numericalFlux = numericalFlux
        self.fluxBoundaryConditions=fluxBoundaryConditions
        self.testUpdateIntegrals=False
        self.testUpdateJacobianIntegrals=False
        self.testCalculateShape=False
        self.testNumericalFlux=False
        self.testNumericalFluxJacobian=False
        #determine if the stabilization term is nonlinear from looking at the coefficient flags
        self.stabilizationIsNonlinear = not ((self.coefficients.mass == None or
                                              self.coefficients.mass == 'linear') and
                                             (self.coefficients.advection == None or
                                              self.coefficients.advection == 'linear') and
                                             (self.coefficients.diffusion == None or
                                              self.coefficients.diffusion == 'constant') and
                                             (self.coefficients.potential == None or
                                              self.coefficients.potential == 'linear') and
                                             (self.coefficients.reaction == None or
                                              self.coefficients.reaction == 'linear'))
        #determine if we need element boundary storage
        self.elementBoundaryIntegrals = ((self.conservativeFlux != None) or 
					 (self.numericalFlux != None) or 
					 (self.fluxBoundaryConditions == 'outFlow') or
                                         (self.fluxBoundaryConditions == 'mixedFlow') or
                                         (self.fluxBoundaryConditions == 'setFlow'))
        self.elementBoundaryIntegrals = True
        #
        #calculate some dimensions
        #
        self.nSpace_global = u.femSpace.referenceFiniteElement.referenceElement.dim
        #we'll use the max dof per element in anticipation of p-refinement
        #this should probably be calculated by the FiniteElementSpace object
        self.nDOF_element = self.u.femSpace.max_nDOF_element
        self.nDOF_element_phi = self.phi.femSpace.max_nDOF_element
        self.nFreeDOF_global = self.dirichletConditions.nFreeDOF_global
        #
        NonlinearEquation.__init__(self,self.nFreeDOF_global)
        #
        #The elementQuadrature dictionary argument allows different quadrature
        #for each integral. The goal of the following blocks is
        #to build a single array of points and an array of weights
        #for each integral the matches these points--zero weight
        #if the given point isn't part of the quadrature rule
        #for that integral.
        #
        #TODO: Decide whether this can be put in QuadTools
        #and whether the phi interpolation points can also
        #be rolled into this
        #
        #First calculate the union of all element quadrature points.
        #
        if self.stabilization != None:
            self.integralKeys += ['stab']
        if self.shockCapturing != None:
            self.integralKeys += ['numDiff']
        quadraturePointSet = set()
        for I in self.integralKeys:
            quadraturePointSet |= set(
                [(p[X],p[Y],p[Z]) for p in elementQuadrature[I].points])
        self.nQuadraturePoints_element = len(quadraturePointSet)
        self.nQuadraturePoints_global = (self.mesh.nElements_global*
                                         self.nQuadraturePoints_element)
        #
        #Now build a dictionary at each element quadrature point which
        #contains the weights for each integral
        #
        # e.g. quadratureWeightDict[p]['m'] is the weight at p for the
        # mass integral
        #
        #initialize all weights to 0
        quadratureWeightDict={}
        for k,p in enumerate(quadraturePointSet):
            quadratureWeightDict[p]={}
            for I in self.integralKeys:
                quadratureWeightDict[p][I]=0.0
        #set the nonzero weights
        for I in self.integralKeys:
            for w,p in zip(elementQuadrature[I].weights,elementQuadrature[I].points):
                quadratureWeightDict[(p[X],p[Y],p[Z])][I]=w
        #
        # Now create the desired point and weight arrays
        #
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
        #storage dictionaries
        #TODO: Cull unused variables
        #
        #element quantities nElements_global x dim object
        #
        self.e={}
        self.scalars_element = ['diameter']
        if self.conservativeFlux != None:
            self.scalars_element += ['conservationResidual']
            if self.conservativeFlux == 'pwc':
                self.scalars_element += ['conservationCorrectionPWC']
        for k in self.scalars_element:
            self.e[k] = Numeric.zeros(
                (self.mesh.nElements_global,),Numeric.Float)
        for eN in range(self.mesh.nElements_global):
            self.e['diameter'][eN] = self.mesh.elementDiametersArray[eN]
        #
        #phi interpolation point quantities nElements_global x nPhiInterpolationPoints_element
        #
        #These are the points required to project phi(u) into the finite element space for phi
        #
        if self.coefficients.potential == 'nonlinear':
            self.phi_ip={}
            self.phi_ip['x']=self.phi.femSpace.interpolationPoints
            self.nPhiInterpolationPoints_element = self.phi.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints
            self.scalars_phi_ip = ['u',
                                   'phi','dphi',
                                   'm','dm',
                                   'r','dr']
            for k in self.scalars_phi_ip: self.phi_ip[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.nPhiInterpolationPoints_element),
                Numeric.Float)
            self.vectors_phi_ip = ['f','df']
            for k in self.vectors_phi_ip: self.phi_ip[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.nPhiInterpolationPoints_element,
                 self.nSpace_global),
                Numeric.Float)
            self.tensors_phi_ip = ['a','da']
            for k in self.tensors_phi_ip: self.phi_ip[k]=Numeric.zeros(
                (self.mesh.nElements_global,
                 self.nPhiInterpolationPoints_element,
                 self.nSpace_global,
                 self.nSpace_global),Numeric.Float)
        #time stepping 
        self.T=0.0
        #
        #quadrature point quantities nElements_global x nQuadraturePoints_element x dim object
        #
        self.q={}
        self.scalars_quadrature = ['u',
                                   'phi','dphi',
                                   'm','dm','mt','dmt','dx_m',
                                   'dx_a','div(a)',
                                   'dx_f','div(f)',
                                   'r','dr','dx_r',
                                   'mt*dx_m',
                                   'dmt*dx_m',
                                   'r*dx_r',
                                   'dr*dx_r',
                                   'det(J)',
                                   'cfl',
                                   'pe']
        if self.stabilization != None:
            self.scalars_quadrature += ['tau','dtau','dx_stab','tau*dx_stab']
        if self.shockCapturing != None:
            self.scalars_quadrature += ['numDiff','dnumDiff','dx_numDiff','numDiff*dx_numDiff']
        if self.stabilization != None or  self.shockCapturing != None:
            self.scalars_quadrature += ['pdeResidual']
        for k in self.scalars_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,self.nQuadraturePoints_element),
            Numeric.Float)
        self.vectors_quadrature = ['grad(u)',
                                   'grad(phi)',
                                   'f','df',
                                   'f*dx_f','df*dx_f',
                                   'velocity']
        for k in self.vectors_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),
            Numeric.Float)
        self.q['x'] = Numeric.zeros((self.mesh.nElements_global,
                                     self.nQuadraturePoints_element,
                                     3),
                                    Numeric.Float)
        self.tensors_quadrature = ['a',
                                   'da',
                                   'a*dx_a',
                                   'da*dx_a',
                                   'J',
                                   'inverse(J)']
        for k in self.tensors_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global,
             self.nSpace_global),
            Numeric.Float)
        self.shape_quadrature = ['v','w']
        if self.stabilization != None or  self.shockCapturing != None:
            self.shape_quadrature += ['dpdeResidual','Lstar*w','ddiv(f)','ddiv(a)']
        for k in self.shape_quadrature: self.q[k]=Numeric.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nDOF_element),Numeric.Float)
        self.gradient_shapeGradient_quadrature = ['grad(u)*grad(w)','grad(phi)*grad(w)']
        for k in self.gradient_shapeGradient_quadrature:
            self.q[k]=Numeric.zeros((self.mesh.nElements_global,
                                     self.nQuadraturePoints_element,
                                     self.nDOF_element,
                                     self.nSpace_global,
                                     self.nSpace_global),Numeric.Float)
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
            self.eb_global = {}
            self.scalars_elementBoundary_global=['fluxAverage','conservationFluxPWC','conservationFluxPWL','heb']
            for k in self.scalars_elementBoundary_global:
                self.eb_global[k] = Numeric.zeros((self.mesh.nElementBoundaries_global,),Numeric.Float)
            if self.conservativeFlux == 'pwl':
                self.eb_global['wFlux'] = Numeric.zeros((self.mesh.nElementBoundaries_global,self.nDOF_element),Numeric.Float)
                self.eb_global['fluxAverage_wFlux'] = Numeric.zeros((self.mesh.nElementBoundaries_global,self.nDOF_element),Numeric.Float)
            self.scalars_elementBoundaryQuadrature_global = ['massAverage','velocityJump','fluxJump','fluxAverage','dx_f','dx_a','advectiveFlux*dx_f','diffusiveFlux*dx_a','penalty']
            for k in self.scalars_elementBoundaryQuadrature_global: self.ebq_global[k]=Numeric.zeros(
                (self.mesh.nElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary),
                Numeric.Float)
            self.vectors_elementBoundaryQuadrature_global = ['massJump','velocityAverage','n']
            if self.conservativeFlux == 'pwc':
                self.vectors_elementBoundaryQuadrature_global += ['conservationVelocityPWC']
            elif self.conservativeFlux == 'pwl':
                self.vectors_elementBoundaryQuadrature_global += ['conservationVelocityPWL']
            for k in self.vectors_elementBoundaryQuadrature_global: self.ebq_global[k]=Numeric.zeros(
                (self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),
                Numeric.Float)
            self.ebq_global['x'] = Numeric.zeros(
                (self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),
                Numeric.Float)
            self.fluxJacobian = Numeric.zeros(
                (self.mesh.nElementBoundaries_global,
                 2,
                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                 self.nDOF_element),
                Numeric.Float)
            #
            #element boundary quadrature point quantities (unique to each element)-- nElements_global x nElementBoundaries_element x nElementBoundaryQuadraturePoints_elementBoundary x dim object
            #
            self.scalars_elementBoundaryQuadrature = ['u',
                                                      'phi','dphi',
                                                      'm','dm',
                                                      'r','dr',
                                                      'sqrt(det(g))']
            for k in self.scalars_elementBoundaryQuadrature: self.ebq[k]=Numeric.zeros(
                (self.mesh.nElements_global, self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary),
                Numeric.Float)
            self.vectors_elementBoundaryQuadrature = ['n',
                                                      'grad(u)',
                                                      'f','df',
                                                      'f*dx_f','df*dx_f',
                                                      'grad(phi)',
                                                      'velocity']
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
            self.tensors_elementBoundaryQuadrature = ['a',
                                                      'da',
                                                      'a*dx_a',
                                                      'da*dx_a',
                                                      'inverse(J)']
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
            self.gradient_shapeGradient_normal_elementBoundaryQuadrature = ['grad(u)*n',
                                                                            'grad(phi)*n']
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
        if self.testNumericalFlux and self.elementBoundaryIntegrals:
            self.fluxtemp = Numeric.array(self.ebq_global['advectiveFlux*dx_f'])
        if self.testNumericalFluxJacobian and self.elementBoundaryIntegrals:
            self.fluxJacobiantemp = Numeric.array(self.fluxJacobian)
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
        #
        # set up solvers and vectors for conservative velocity field calculation
        #
        if self.conservativeFlux == 'pwc':
            print "Allocating pwc conservative flux solver"
            self.updateConservationJacobian = True
            self.conservationJacobianPWC = Mat(self.mesh.nElements_global,
                                                     self.mesh.nElements_global,
                                                     self.mesh.nElements_global)
            self.conservationSolver = DenseLU(self.conservationJacobianPWC)
        elif self.conservativeFlux == 'pwl':
            print "Allocating pwl conservative flux solver and data structures"
            self.updateConservationJacobian = True
            self.globalNode2globalElementList = self.mesh.nodeElementsList
            self.globalNodeGlobalElement2StarElement = []
            self.subdomainL=[]
            self.subdomainV=[]
            self.subdomainR=[]
            self.conservationSolverPWL=[]
	    self.nodeStarElementsArray = Numeric.zeros((self.mesh.nElements_global,
							self.mesh.nNodes_element),Numeric.Int)
	    self.nodeStarElementNeighborsArray = Numeric.zeros((self.mesh.nElements_global,
								self.mesh.nNodes_element,
								self.mesh.nElementBoundaries_element),
							       Numeric.Int)
            self.nodeStarOffset = Numeric.zeros((self.mesh.nNodes_global,),Numeric.Int)
            self.nodeStarJacobianOffset = Numeric.zeros((self.mesh.nNodes_global,),Numeric.Int)
            offset=0
            jacobianOffset=0
            for I in range(self.mesh.nNodes_global):
                self.globalNode2globalElementList[I].sort()
                self.globalNodeGlobalElement2StarElement.append(dict([(eN_global,eN_node) for eN_node,eN_global in enumerate(self.globalNode2globalElementList[I])]))
                self.nodeStarOffset[I]=offset
                self.nodeStarJacobianOffset[I] = jacobianOffset
                offset += self.mesh.nElements_node[I]
                jacobianOffset += self.mesh.nElements_node[I]**2
                self.subdomainR.append(Vec(self.mesh.nElements_node[I]))
                self.subdomainV.append(Vec(self.mesh.nElements_node[I]))
                self.subdomainL.append(Mat(self.mesh.nElements_node[I],self.mesh.nElements_node[I]))
                self.conservationSolverPWL.append(LinearSolvers.DenseLU(self.subdomainL[I]))
	    self.starU = Numeric.zeros((offset,),Numeric.Float)
	    self.starR = Numeric.zeros((offset,),Numeric.Float)
            self.starJacobian = Numeric.zeros((jacobianOffset,),Numeric.Float)
	    for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
		ebN = self.mesh.interiorElementBoundariesArray[ebNI]
		left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
		right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
		left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
		right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                for i in range(self.mesh.nNodes_element):
                    left_I = self.mesh.elementNodesArray[left_eN_global,i]
		    self.nodeStarElementsArray[left_eN_global,i] = self.globalNodeGlobalElement2StarElement[left_I][left_eN_global]
                    if i != left_ebN_element:
			self.nodeStarElementNeighborsArray[left_eN_global,i,left_ebN_element] = self.globalNodeGlobalElement2StarElement[left_I][right_eN_global]
	    for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
		ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
		eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
		ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
		for i in range(self.mesh.nNodes_element):
		    I = self.mesh.elementNodesArray[eN_global,i]
		    self.nodeStarElementsArray[eN_global,i] = self.globalNodeGlobalElement2StarElement[I][eN_global]
        self.inflowBoundary = Numeric.zeros((self.mesh.nExteriorElementBoundaries_global,),Numeric.Int)
        self.inflowBoundaryBC = Numeric.zeros((self.mesh.nExteriorElementBoundaries_global,),Numeric.Int)
        self.inflowBoundaryBC_values = Numeric.zeros((self.mesh.nExteriorElementBoundaries_global,self.nDOF_element),Numeric.Float)
        self.inflowFlux = Numeric.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),Numeric.Float)
        self.internalNodes = set(range(self.mesh.nNodes_global))
	#tag internal nodes this is ought to be in mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global,i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray=Numeric.zeros((self.nNodes_internal,),Numeric.Int)
        for nI,n in enumerate(self.internalNodes):
            self.internalNodesArray[nI]=n
        print "updating local to global mapping matrix"
        self.updateLocal2Global()
        print "building time inegration object"
        print TimeIntegrationClass
        #mwf added extra arguments?
        if TimeIntegrationClass == SSPRKintegration:
            self.timeIntegration = TimeIntegrationClass(self,tOrder)
        else:
            self.timeIntegration = TimeIntegrationClass(self)
        #end if
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
        #mwf debug
        #print 'after getResidual elementRes= \n',self.elementResidual
        #print 'residual= \n',r
        #print "element residual \n" + `self.elementResidual` + "\n r \n"+`r`
        return r
    def getJacobian(self,jacobian):
	if self.matType == 'csr':
	    femIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
					  jacobian)
	else:
	    jacobian.flat[:]=0.0
	self.updateElementJacobian()
	self.updateElementBoundaryJacobian()
	if self.matType == 'csr':
	    self.getJacobian_CSR(jacobian)
	else:
	    self.getJacobian_dense(jacobian)
        #print "element jacobian \n" + `self.elementJacobian` + "\n jacobian \n"+`jacobian`
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
        if self.numericalFlux != None:
            femIntegrals.updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(self.mesh.nInteriorElementBoundaries_global,
											   self.nDOF_element,
											   self.nElementBoundaryQuadraturePoints_elementBoundary,
											   self.mesh.nElementBoundaries_element,
											   self.nFreeDOF_global,
											   self.mesh.interiorElementBoundariesArray,
											   self.mesh.elementBoundaryElementsArray,
											   self.mesh.elementBoundaryLocalElementBoundariesArray,
											   self.l2g['nFreeDOF'],
											   self.l2g['freeLocal'],
											   self.l2g['freeGlobal'],
											   self.fluxJacobian,
											   self.ebq['w'],
											   jacobian)
            femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(self.mesh.nExteriorElementBoundaries_global,
											   self.nDOF_element,
											   self.nElementBoundaryQuadraturePoints_elementBoundary,
											   self.mesh.nElementBoundaries_element,
											   self.nFreeDOF_global,
											   self.mesh.exteriorElementBoundariesArray,
											   self.mesh.elementBoundaryElementsArray,
											   self.mesh.elementBoundaryLocalElementBoundariesArray,
											   self.l2g['nFreeDOF'],
											   self.l2g['freeLocal'],
											   self.l2g['freeGlobal'],
											   self.fluxJacobian,
											   self.ebq['w'],
											   jacobian)
        elif ( self.fluxBoundaryConditions == 'outFlow' or self.fluxBoundaryConditions == 'mixedFlow'):
	    femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(self.mesh.nExteriorElementBoundaries_global,
											   self.nDOF_element,
											   self.nElementBoundaryQuadraturePoints_elementBoundary,
											   self.mesh.nElementBoundaries_element,
											   self.nFreeDOF_global,
											   self.mesh.exteriorElementBoundariesArray,
											   self.mesh.elementBoundaryElementsArray,
											   self.mesh.elementBoundaryLocalElementBoundariesArray,
											   self.l2g['nFreeDOF'],
											   self.l2g['freeLocal'],
											   self.l2g['freeGlobal'],
											   self.fluxJacobian,
											   self.ebq['w'],
											   jacobian)
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
        #self.initializeJacobian()
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
        if self.numericalFlux != None:
            femIntegrals.updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(self.mesh.nInteriorElementBoundaries_global,
                                                                                         self.nDOF_element,
                                                                                         self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                         self.mesh.nElementBoundaries_element,
                                                                                         self.mesh.interiorElementBoundariesArray,
                                                                                         self.mesh.elementBoundaryElementsArray,
                                                                                         self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                         self.l2g['nFreeDOF'],
                                                                                         self.l2g['freeLocal'],
                                                                                         self.csrRowIndeces,
                                                                                         self.csrColumnOffsets_eb,
                                                                                         self.fluxJacobian,
											 self.ebq['w'],
                                                                                         jacobian)
            femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(self.mesh.nExteriorElementBoundaries_global,
                                                                                         self.nDOF_element,
                                                                                         self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                         self.mesh.nElementBoundaries_element,
                                                                                         self.mesh.exteriorElementBoundariesArray,
                                                                                         self.mesh.elementBoundaryElementsArray,
                                                                                         self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                         self.l2g['nFreeDOF'],
                                                                                         self.l2g['freeLocal'],
                                                                                         self.csrRowIndeces,
                                                                                         self.csrColumnOffsets_eb,
                                                                                         self.fluxJacobian,
                                                                                         self.ebq['w'],
                                                                                         jacobian)
        elif (self.fluxBoundaryConditions == 'outFlow' or self.fluxBoundaryConditions == 'mixedFlow'):
           femIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(self.mesh.nExteriorElementBoundaries_global,
                                                                                        self.nDOF_element,
                                                                                        self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                        self.mesh.nElementBoundaries_element,
                                                                                        self.mesh.exteriorElementBoundariesArray,
                                                                                        self.mesh.elementBoundaryElementsArray,
                                                                                        self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                        self.l2g['nFreeDOF'],
                                                                                        self.l2g['freeLocal'],
                                                                                        self.csrRowIndeces,
                                                                                        self.csrColumnOffsets_eb,
                                                                                        self.fluxJacobian,
                                                                                        self.ebq['w'],
                                                                                        jacobian)
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
        return jacobian
    def updateElementResidual(self):
        """Calculate all the element residuals"""
        self.elementResidual.flat[:]=0.0
        #mwf debug
        #print 'entering elemRes= ',self.elementResidual
        #ara optimization idea
#          femIntegrals.updateMassAdvectionDiffusionReaction(self.mesh.nElements_global,
#                                   self.nQuadraturePoints_element,
#                                   self.nDOF_element,
#                                   self.nSpace_global,
#                                   self.coefficients.mass != None,
#                                   self.coefficients.advection != None,
#                                   self.coefficients.diffusion != None,
#                                   self.coefficients.reaction != None,
#                                   self.q['mt*dx_m'],
#                                   self.q['w'],
#                                   self.q['f*dx_f'],
#                                   self.q['grad(w)'],
#                                   self.q['a*dx_a'],
#                                   self.q['grad(phi)*grad(w)'],
#                                   self.q['r*dx_r'],
#                                   self.elementResidual)
        if self.coefficients.mass != None:
            self.updateMass(self.elementResidual)
            #mwf debug
            #print 'after updateMass elemRes= ',self.elementResidual
        if self.coefficients.advection != None:
            self.updateAdvection(self.elementResidual)
            #mwf debug
            #print 'after updateAdvection elemRes= ',self.elementResidual
        if self.coefficients.diffusion != None:
            self.updateDiffusion(self.elementResidual)
            #mwf debug
            #print 'after updateDiffusion elemRes= ',self.elementResidual
        if self.coefficients.reaction != None:
            self.updateReaction(self.elementResidual)
            #mwf debug
            #print 'after updateReaction elemRes= ',self.elementResidual
        if self.stabilization != None:
            self.updateStabilization(self.elementResidual)
            #mwf debug
            #print 'after updateStabilization elemRes= ',self.elementResidual
        if self.shockCapturing != None:
            self.updateShockCapturing(self.elementResidual)
            #mwf debug
            #print 'after updateShockCapturing elemRes= ',self.elementResidual
        if self.numericalFlux != None:
            if self.coefficients.advection != None:
                self.updateAdvectionInteriorElementBoundaryFlux(self.elementResidual)
            if self.coefficients.diffusion != None:
                self.updateDiffusionInteriorElementBoundaryFlux(self.elementResidual)
            #mwf debug
            #print 'after updateBoundaryCoeff elemRes= ',self.elementResidual
        if self.fluxBoundaryConditions == 'outFlow':
            if self.coefficients.advection != None:
                self.updateAdvectionExteriorElementBoundaryFlux(self.elementResidual)
        if (self.fluxBoundaryConditions == 'mixedFlow' or
            self.fluxBoundaryConditions == 'setFlow'):
            if self.coefficients.advection != None:
                self.updateAdvectionExteriorElementBoundaryFlux(self.elementResidual)
            if self.coefficients.diffusion != None:
                self.updateDiffusionExteriorElementBoundaryFlux(self.elementResidual)
        if self.dirichletNodeSetList != None:
            for eN in range(self.mesh.nElements_global):
                for j in self.dirichletNodeSetList[eN]:
                    J = self.trialSpace.dofMap.l2g[eN,j]
                    self.u.dof[J] = self.dirichletValues[(eN,j)]
                    self.elementResidual[eN,j] = self.u.dof[J]-self.dirichletValues[(eN,j)]
                    
    def updateElementJacobian(self):
        self.elementJacobian.flat[:]=0.0
#          femIntegrals.updateMassAdvectionDiffusionReactionJacobian(self.mesh.nElements_global,
#                                          self.nQuadraturePoints_element,
#                                          self.nDOF_element,
#                                          self.nSpace_global,
#                                          self.coefficients.mass != None and self.timeIntegration.massIsImplicit,
#                                          self.coefficients.advection != None and self.timeIntegration.advectionIsImplicit,
#                                          self.coefficients.diffusion != None and self.timeIntegration.diffusionIsImplicit,
#                                          self.coefficients.reaction != None and self.timeIntegration.reactionIsImplicit,
#                                          self.q['dmt*dx_m'],
#                                          self.q['v*w'],
#                                          self.q['df*dx_f'],
#                                          self.q['v*grad(w)'],
#                                          self.phi.femSpace.dofMap.l2g,
#                                          self.q['a*dx_a'],
#                                          self.q['da*dx_a'],
#                                          self.q['grad(phi)*grad(w)'],
#                                          self.dphi.dof,
#                                          self.q['v'],
#                                          self.q['grad(v)*grad(w)'],
#                                          self.q['dr*dx_r'],
#                                          self.q['v*w'],
#                                          self.elementJacobian)
        if self.coefficients.mass != None:
            if self.timeIntegration.massIsImplicit:
                self.updateMassJacobian(self.elementJacobian)
        if self.coefficients.advection != None:
            if self.timeIntegration.advectionIsImplicit:
                self.updateAdvectionJacobian(self.elementJacobian)
        if self.coefficients.diffusion != None:
            if self.timeIntegration.diffusionIsImplicit:
                self.updateDiffusionJacobian(self.elementJacobian)
        if self.coefficients.reaction != None:
            if self.timeIntegration.reactionIsImplicit:
                self.updateReactionJacobian(self.elementJacobian)
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
        #Note: The element boundary fluxes yield global jacobian entries
        #so they're in getJacobian
    def updateElementBoundaryJacobian(self):
        if self.numericalFlux != None or self.fluxBoundaryConditions == 'outFlow':
            self.fluxJacobian.flat[:]=0.0
        if self.numericalFlux != None:
            if self.coefficients.advection != None:
                self.updateInteriorNumericalAdvectiveFluxJacobian(self.fluxJacobian)
            if self.coefficients.diffusion != None:
                self.updateInteriorNumericalDiffusiveFluxJacobian(self.fluxJacobian)
        if self.fluxBoundaryConditions == 'outFlow':
            if self.coefficients.advection != None:
                self.updateExteriorNumericalAdvectiveFluxJacobian(self.fluxJacobian)
        if self.fluxBoundaryConditions == 'mixedFlow':
            if self.coefficients.advection != None:
                self.updateExteriorNumericalAdvectiveFluxJacobian(self.fluxJacobian)
            if self.coefficients.advection != None:
                self.updateExteriorNumericalDiffusiveFluxJacobian(self.fluxJacobian)
        if self.numericalFlux != None or self.fluxBoundaryConditions == 'outFlow' or self.fluxBoundaryConditions == 'mixedFlow':
            self.timeIntegration.updateNumericalFluxJacobian(self.fluxJacobian)
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
        #
        #get functions of (t,x,u) at the quadrature points
        #
        self.coefficients.evaluate(t    = self.T,
                                   x    = Numeric.reshape(
            self.q['x'].flat,
            (self.nQuadraturePoints_global,3)),
                                   u    = self.q['u'].flat,
                                   m    = self.q['m'].flat,
                                   dm   = self.q['dm'].flat,
                                   f    = Numeric.reshape(
            self.q['f'].flat,(self.nQuadraturePoints_global,self.nSpace_global)),
                                   df   = Numeric.reshape(
            self.q['df'].flat,(self.nQuadraturePoints_global,self.nSpace_global)),
                                   a    = Numeric.reshape(
            self.q['a'].flat,(self.nQuadraturePoints_global,self.nSpace_global,self.nSpace_global)),
                                   da   = Numeric.reshape(
            self.q['da'].flat,(self.nQuadraturePoints_global,self.nSpace_global,self.nSpace_global)),
                                   phi  = self.q['phi'].flat,
                                   dphi = self.q['dphi'].flat,
                                   r    = self.q['r'].flat,
                                   dr   = self.q['dr'].flat)
        #print "q \n"+`self.q`
        #
        # m_t at the quadrature points
        #
        self.timeIntegration.updateMass(self.q['m'],self.q['mt'],self.q['dm'],self.q['dmt'])
        #mwf debug
        #print 'ScalarTransport out of timeInt.updateMass '
        #print '\t q[m]=\n ',self.q['m']
        #print '\t q[mt]=\n ',self.q['mt']
        #print '\t q[dm]=\n ',self.q['dm']
        #print '\t q[dmt]=\n ',self.q['dmt']
        #print '\t q[f]= \n',self.q['f']
        #print '\t q[df]= \n',self.q['df']
        #
        #
        #calculate div(f),div(a),tau, and Lstar for use in residual and stabilization calculations
        #
        if self.stabilization != None or self.shockCapturing != None:
            femIntegrals.calculateDiv_f(self.mesh.nElements_global,
                                        self.nQuadraturePoints_element,
                                        self.nDOF_element,
                                        self.nSpace_global,
                                        self.q['df'],
                                        self.q['grad(u)'],
                                        self.q['grad(v)'],
                                        self.q['div(f)'],
                                        self.q['ddiv(f)']);
            femIntegrals.calculateDiv_a(self.mesh.nElements_global,
                                        self.nQuadraturePoints_element,
                                        self.nDOF_element,
                                        self.nSpace_global,
                                        self.phi.femSpace.dofMap.l2g,
                                        self.q['da'],
                                        self.dphi.dof,
                                        self.q['grad(phi)'],
                                        self.q['grad(u)'],
                                        self.q['grad(v)'],
                                        self.q['div(a)'],
                                        self.q['ddiv(a)']);
        if self.stabilization != None:
            femIntegrals.calculateStabilizationADR(self.mesh.nElements_global,
                                                self.nQuadraturePoints_element,
                                                self.nSpace_global,
                                                self.stabilization,
                                                self.e['diameter'],
                                                self.q['df'],
                                                self.q['a'],
                                                self.q['da'],
                                                self.q['grad(phi)'],
                                                self.q['dphi'],
                                                self.q['dr'],
                                                self.q['dmt'],
                                                self.q['pe'],
                                                self.q['cfl'],
                                                self.q['tau'])
            #
            # L^* w at the quadrature points
            #
            self.calculateAdjoint()
        else:
            femIntegrals.calculateDimensionlessNumbersADR(self.mesh.nElements_global,
                                                          self.nQuadraturePoints_element,
                                                          self.nSpace_global,
                                                          self.e['diameter'],
                                                          self.q['df'],
                                                          self.q['a'],
                                                          self.q['dphi'],
                                                          self.q['dr'],
                                                          self.q['dmt'],
                                                          self.q['pe'],
                                                          self.q['cfl'])
        #
        # get the rest of the terms at the quadrature points allowing time integration
        # to modify
        #
        if self.stabilization != None:
            self.timeIntegration.updateStabilization(self.q['tau'])
            self.timeIntegration.updateAdjoint(self.q['Lstar*w'])
        #
        #phi
        #
        #
        #get functions of (t,x,u) at phi_ip if necessary
        #
        if (self.coefficients.potential  == 'nonlinear' and self.coefficients.diffusion != None):
            #the commented out line is what we should use in general to calculation u at the interpolation points
#          femIntegrals.calculateFiniteElementFunction(self.mesh.nElements_global,
#                                                      self.nPhiInterpolationPoints_element,
#                                                      self.nDOF_element,
#                                                      self.nSpace_global,self.u.femSpace.dofMap.l2g,self.u.dof,
#                                                      self.phi_ip['v'],self.phi_ip['grad(v)'],self.phi_ip['grad(v)*grad(w)'],
#                                                      self.phi_ip['u'],self.phi_ip['grad(u)'],self.phi_ip['grad(u)*grad(w)'])
            #for now use the fact that the phi_ip/dof are the u_ip/dof
            for eN in range(self.mesh.nElements_global):
                for j in range(self.nDOF_element):
                    J = self.u.femSpace.dofMap.l2g[eN,j]
                    self.phi_ip['u'][eN,j] = self.u.dof[J]
            self.coefficients.evaluate(t    = self.T,
                                       x    = Numeric.reshape(self.phi_ip['x'].flat,
                                                      (self.mesh.nElements_global*self.nPhiInterpolationPoints_element,self.nSpace_global)),
                                       u    = self.phi_ip['u'].flat,
                                       m    = self.phi_ip['m'].flat,
                                       dm   = self.phi_ip['dm'].flat,
                                       f    = Numeric.reshape(self.phi_ip['f'].flat,
                                                      (self.mesh.nElements_global*self.nPhiInterpolationPoints_element,self.nSpace_global)),
                                       df   = Numeric.reshape(self.phi_ip['df'].flat,
                                                      (self.mesh.nElements_global*self.nPhiInterpolationPoints_element,self.nSpace_global)),
                                       a    = Numeric.reshape(self.phi_ip['a'].flat,
                                                      (self.mesh.nElements_global*self.nPhiInterpolationPoints_element,self.nSpace_global,self.nSpace_global)),
                                       da   = Numeric.reshape(self.phi_ip['da'].flat,
                                                      (self.mesh.nElements_global*self.nPhiInterpolationPoints_element,self.nSpace_global,self.nSpace_global)),
                                       phi  = self.phi_ip['phi'].flat,
                                       dphi = self.phi_ip['dphi'].flat,
                                       r    = self.phi_ip['r'].flat,
                                       dr   = self.phi_ip['dr'].flat)
            self.phi.projectFromInterpolationConditions(self.phi_ip['phi'])
            self.dphi.projectFromInterpolationConditions(self.phi_ip['dphi'])
            self.calculateFiniteElementPhi()
        else:
            self.phi.dof.flat[:]=self.u.dof.flat
            self.dphi.dof.flat[:]=1.0
            self.q['phi'].flat[:]=self.q['u'].flat
            self.q['dphi'].flat[:]=1.0
            self.q['grad(phi)'].flat[:]=self.q['grad(u)'].flat
            self.q['grad(phi)*grad(w)'].flat[:]=self.q['grad(u)*grad(w)'].flat
        self.timeIntegration.updateAdvection(self.q['f'],self.q['df'])
        self.timeIntegration.updateDiffusion(self.q['a'],self.q['da'])
        self.timeIntegration.updateGradients(self.dphi.dof,self.q['grad(phi)*grad(w)'],self.q['grad(u)*grad(w)'])
        self.timeIntegration.updateReaction(self.q['r'],self.q['dr'])
        if self.stabilization != None or self.shockCapturing != None:
            self.timeIntegration.updateDivF(self.q['div(f)'],self.q['ddiv(f)'])
            self.timeIntegration.updateDivA(self.q['div(a)'],self.q['ddiv(a)'])
            femIntegrals.calculatePDEResidualADR(self.mesh.nElements_global,
                                                 self.nQuadraturePoints_element,
                                                 self.q['div(f)'],
                                                 self.q['div(a)'],
                                                 self.q['r'],
                                                 self.q['mt'],
                                                 self.q['pdeResidual'])
            femIntegrals.calculatePDEResidualJacobianADR(self.mesh.nElements_global,
                                                         self.nQuadraturePoints_element,
                                                         self.nDOF_element,
                                                         self.q['ddiv(f)'],
                                                         self.q['ddiv(a)'],
                                                         self.q['dr'],
                                                         self.q['dmt'],
                                                         self.q['v'],
                                                         self.q['dpdeResidual'])
        if self.shockCapturing != None:
            femIntegrals.calculateShockCapturingADR(self.mesh.nElements_global,
                                                    self.nQuadraturePoints_element,
                                                    self.nSpace_global,
                                                    self.shockCapturing,
                                                    self.shockCapturingDiffusion,
                                                    self.e['diameter'],
                                                    self.q['pdeResidual'],
                                                    self.q['mt'],
                                                    self.q['grad(phi)'],
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
        if self.coefficients.diffusion != None:
            femIntegrals.calculateTensorScalarProduct(self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.nSpace_global,
                                                      self.q['a'],
                                                      self.q['dx_a'],
                                                      self.q['a*dx_a'])
            if self.timeIntegration.diffusionIsImplicit:
                femIntegrals.calculateTensorScalarProduct(self.mesh.nElements_global,
                                                          self.nQuadraturePoints_element,
                                                          self.nSpace_global,
                                                          self.q['da'],
                                                          self.q['dx_a'],
                                                          self.q['da*dx_a'])
        if self.coefficients.advection != None:
            femIntegrals.calculateVectorScalarProduct(self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.nSpace_global,
                                                      self.q['f'],
                                                      self.q['dx_f'],
                                                      self.q['f*dx_f'])
            if self.timeIntegration.advectionIsImplicit:
                femIntegrals.calculateVectorScalarProduct(self.mesh.nElements_global,
                                                          self.nQuadraturePoints_element,
                                                          self.nSpace_global,
                                                          self.q['df'],
                                                          self.q['dx_f'],
                                                          self.q['df*dx_f'])
        if self.coefficients.reaction != None:
            femIntegrals.calculateScalarScalarProduct(self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.q['r'],
                                                      self.q['dx_r'],
                                                      self.q['r*dx_r'])
            if self.timeIntegration.reactionIsImplicit:
                femIntegrals.calculateScalarScalarProduct(self.mesh.nElements_global,
                                                          self.nQuadraturePoints_element,
                                                          self.q['dr'],
                                                          self.q['dx_r'],
                                                          self.q['dr*dx_r'])
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
        #
        #get coefficients at the element boundary quadrature points
        #
        self.coefficients.evaluate(t    = self.T,
                                   x    = Numeric.reshape(
            self.ebq['x'].flat,
            (self.nElementBoundaryQuadraturePoints_global,3)),
                                   u    = self.ebq['u'].flat,
                                   m    = self.ebq['m'].flat,
                                   dm   = self.ebq['dm'].flat,
                                   f    = Numeric.reshape(
            self.ebq['f'].flat,(self.nElementBoundaryQuadraturePoints_global,self.nSpace_global)),
                                   df   = Numeric.reshape(
            self.ebq['df'].flat,(self.nElementBoundaryQuadraturePoints_global,self.nSpace_global)),
                                   a    = Numeric.reshape(
            self.ebq['a'].flat,(self.nElementBoundaryQuadraturePoints_global,self.nSpace_global,self.nSpace_global)),
                                   da   = Numeric.reshape(
            self.ebq['da'].flat,(self.nElementBoundaryQuadraturePoints_global,self.nSpace_global,self.nSpace_global)),
                                   phi  = self.ebq['phi'].flat,
                                   dphi = self.ebq['dphi'].flat,
                                   r    = self.ebq['r'].flat,
                                   dr   = self.ebq['dr'].flat)
        #
        #get phi at the quadrature points if necessary
        #
        if self.coefficients.potential  == 'nonlinear':
            self.calculateFiniteElementPhi_elementBoundaries()
        else:
            self.ebq['phi'].flat[:] =self.ebq['u'].flat
            self.ebq['grad(phi)'].flat[:] = self.ebq['grad(u)'].flat
        #
        # calculate the averages and jumps at element boundaries
        #
        # TODO: make sure we need all this info
        #
        #if self.conservativeFlux != None:
        #    femIntegrals.calculateInteriorElementBoundaryVelocityAverage()
        #if self.conservativeFlux != None:
        #    femIntegrals.calculateExteriorElementBoundaryVelocityAverage()
        if self.numericalFlux != None:
            if self.coefficients.advection != None:
                self.calculateInteriorNumericalAdvectiveFlux()
            if self.coefficients.diffusion != None:
                self.calculateInteriorNumericalDiffusiveFlux()        
        if self.fluxBoundaryConditions == 'outFlow':
            if self.coefficients.advection != None:
                self.calculateExteriorNumericalAdvectiveFlux()
        if self.fluxBoundaryConditions == 'mixedFlow':
            if self.coefficients.advection != None:
                self.calculateExteriorNumericalAdvectiveFlux()
            if self.coefficients.advection != None:
                self.calculateExteriorNumericalDiffusiveFlux()
        if self.fluxBoundaryConditions == 'setFlow':
            self.flowBC_setter.setExteriorFlux(self.mesh.nExteriorElementBoundaries_global,
                                               self.nElementBoundaryQuadraturePoints_elementBoundary,
                                               self.mesh.nElementBoundaries_element,
                                               self.nSpace_global,
                                               self.mesh.exteriorElementBoundariesArray,
                                               self.mesh.elementBoundaryElementsArray,
                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
                                               self.ebq['x'],
                                               self.ebq['n'],
                                               self.ebq_global['dx_f'],
                                               self.ebq_global['dx_a'],
                                               self.ebq_global['advectiveFlux*dx_f'],
                                               self.ebq_global['diffusiveFlux*dx_a'])            
        self.timeIntegration.updateNumericalFlux(self.ebq_global['advectiveFlux*dx_f'],
                                                 self.ebq_global['diffusiveFlux*dx_a'])
	if self.conservativeFlux != None:
	    femIntegrals.calculateInteriorElementBoundaryVelocities(self.mesh.nInteriorElementBoundaries_global,
								    self.nElementBoundaryQuadraturePoints_elementBoundary,
								    self.mesh.nElementBoundaries_element,
								    self.nSpace_global,
								    self.mesh.interiorElementBoundariesArray,
								    self.mesh.elementBoundaryElementsArray,
								    self.mesh.elementBoundaryLocalElementBoundariesArray,
								    self.ebq['m'],
								    self.ebq['a'],
								    self.ebq['grad(phi)'],
								    self.ebq['f'],
								    self.ebq_global['velocityAverage'],
								    self.ebq_global['velocityJump'],
								    self.ebq_global['massAverage'],
								    self.ebq_global['massJump'])
	    femIntegrals.calculateExteriorElementBoundaryVelocities(self.mesh.nExteriorElementBoundaries_global,
								    self.nElementBoundaryQuadraturePoints_elementBoundary,
								    self.mesh.nElementBoundaries_element,
								    self.nSpace_global,
								    self.mesh.exteriorElementBoundariesArray,
								    self.mesh.elementBoundaryElementsArray,
								    self.mesh.elementBoundaryLocalElementBoundariesArray,
								    self.ebq['m'],
								    self.ebq['a'],
								    self.ebq['grad(phi)'],
								    self.ebq['f'],
								    self.ebq_global['velocityAverage'],
								    self.ebq_global['velocityJump'],
								    self.ebq_global['massAverage'],
								    self.ebq_global['massJump'])
#         self.eb_global['fluxAverage'].flat[:]=0.0
#         for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
#             ebN = self.mesh.interiorElementBoundariesArray[ebNI]
#             left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#             right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
#             left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#             right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
#             left_normals   = self.ebq['n'][left_eN_global,left_ebN_element]
#             right_normals  = self.ebq['n'][right_eN_global,right_ebN_element]
#             for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                 self.ebq_global['velocityAverage'][ebN,k,:] = 0.5*(
#                     self.ebq['f'][left_eN_global,left_ebN_element,k] -
#                     tenvec(self.ebq['a'][left_eN_global,left_ebN_element,k],self.ebq['grad(phi)'][left_eN_global,left_ebN_element,k])
#                     +
#                     self.ebq['f'][right_eN_global,right_ebN_element,k] -
#                     tenvec(self.ebq['a'][right_eN_global,right_ebN_element,k],self.ebq['grad(phi)'][right_eN_global,right_ebN_element,k]))
#                 self.ebq_global['velocityJump'][ebN,k] = (
#                     dot(self.ebq['f'][left_eN_global,left_ebN_element,k] -
#                          tenvec(self.ebq['a'][left_eN_global,left_ebN_element,k],self.ebq['grad(phi)'][left_eN_global,left_ebN_element,k]),
#                          left_normals[k])
#                     +
#                     dot(self.ebq['f'][right_eN_global,right_ebN_element,k] -
#                          tenvec(self.ebq['a'][right_eN_global,right_ebN_element,k],self.ebq['grad(phi)'][right_eN_global,right_ebN_element,k]),
#                          right_normals[k]))
#                 self.ebq_global['massAverage'][ebN,k] = 0.5*(
#                     self.ebq['m'][left_eN_global,left_ebN_element,k]
#                     +
#                     self.ebq['m'][right_eN_global,right_ebN_element,k])
#                 self.ebq_global['massJump'][ebN,k,:] = (
#                     self.ebq['m'][left_eN_global,left_ebN_element,k]*left_normals[k]
#                     +
#                     self.ebq['m'][right_eN_global,right_ebN_element,k]*right_normals[k])
#                 self.eb_global['fluxAverage'][ebN] += dot(self.ebq_global['velocityAverage'][ebN,k],left_normals[k])*self.ebq_global['dx_f'][ebN,k]
#         for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#             ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#             eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#             ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#             normals   = self.ebq['n'][eN_global,ebN_element]
#             for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                 self.ebq_global['velocityAverage'][ebN,k] = (
#                     self.ebq['f'][eN_global,ebN_element,k] -
#                     tenvec(self.ebq['a'][eN_global,ebN_element,k],self.ebq['grad(phi)'][eN_global,ebN_element,k]))
#                 self.ebq_global['velocityJump'][ebN,k] = (
#                     dot(self.ebq['f'][eN_global,ebN_element,k] -
#                          tenvec(self.ebq['a'][eN_global,ebN_element,k],self.ebq['grad(phi)'][eN_global,ebN_element,k]),
#                          normals[k]))
#                 self.ebq_global['massAverage'][ebN,k] = (
#                     self.ebq['m'][eN_global,ebN_element,k])
#                 self.ebq_global['massJump'][ebN,k] = (
#                     self.ebq['m'][eN_global,ebN_element,k]*normals[k])
#                 self.eb_global['fluxAverage'][ebN] += dot(self.ebq_global['velocityAverage'][ebN,k],normals[k])*self.ebq_global['dx_f'][ebN,k]
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
#              print self.q['grad(v)*w']
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
#              for eN in range(self.mesh.nElements_global):
#                  for k,xi in enumerate(self.quadraturePoints):
#                      self.q['dx_'+integral][eN,k] = abs(self.q['det(J)'][eN,k])*self.quadratureWeights[integral][k]
        #
        #get the interpolation points for building phi from phi(u)
        #
        self.phi.femSpace.updateInterpolationPoints()
        self.dphi.femSpace.interpolationPoints = self.phi.femSpace.interpolationPoints
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
        self.u.femSpace.elementMaps.getInverseValuesTrace(self.ebq['inverse(J)'],self.ebq['x'],self.ebq['hat(x)'])
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
        #tensor products of shape information
        for eN in range(self.mesh.nElements_global):
            for ebN in range(self.mesh.nElementBoundaries_element):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    for i in range(self.nDOF_element): 
                        for j in range(self.nDOF_element):
                            self.ebq['v*w'][eN,ebN,k,j,i] = (self.ebq['v'][eN,ebN,k,j]*
                                                                 self.ebq['w'][eN,ebN,k,i])
                        for l in range(self.nSpace_global):
                            #unused storage?
                            self.ebq['grad(v)*w'][eN,ebN,k,j,i,l] = self.ebq['grad(v)'][eN,ebN,k,j,l]*self.ebq['w'][eN,ebN,k,i]
        #quadrature weights for each integral
        for integral in self.elementBoundaryIntegralKeys:
            for ebN in range(self.mesh.nElementBoundaries_global):
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['dx_'+integral][ebN,k] = abs(self.ebq['sqrt(det(g))'][left_eN_global,left_ebN_element,k])*self.elementBoundaryQuadratureWeights[integral][k]
        #element boundary area
        self.eb_global['heb'].flat[:]=0.0
        for ebN in range(self.mesh.nElementBoundaries_global):
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                self.eb_global['heb'][ebN] += self.ebq_global['dx_a'][ebN,k]
                self.ebq_global['penalty'][ebN,k] = 2.0/self.ebq_global['dx_a'][ebN,k]
            self.ebq_global['penalty'][ebN,k] = 2.0/self.eb_global['heb'][ebN]
    def calculateAdjoint(self):
        """
        Calculate the action of the adjoint of the linearized spatial
        operator on the test functions.
        """
        femIntegrals.calculateAdjointADR(self.mesh.nElements_global,
                                         self.nQuadraturePoints_element,
                                         self.nDOF_element,
                                         self.nSpace_global,
                                         self.q['w'],
                                         self.q['grad(w)'],
                                         self.q['df'],
                                         self.q['da'],
                                         self.q['grad(phi)'],
                                         self.q['dr'],
                                         self.q['Lstar*w'])
    def calculateFiniteElementU(self):
        """
        Calculate u, grad(u), and the (outer) product grad(u) x grad(w)
        """
        self.u.getValues(self.q['v'],
                         self.q['u'])
        self.u.getGradientValues(self.q['grad(v)'],
                                 self.q['grad(u)'])
        if  self.coefficients.diffusion != None or self.shockCapturing != None:
            self.u.getGradientTensorValues(self.q['grad(v)*grad(w)'],
                                           self.q['grad(u)*grad(w)'])
        #mwf debug just try to calculate q['u'] manually too
        #manu = Numeric.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),
        #                     Numeric.Float)
        #for eN in range(self.mesh.nElements_global):
        #    for k in range(self.nQuadraturePoints_element):
        #        for j in range(self.nDOF_element):
        #            J = self.u.femSpace.dofMap.l2g[eN,j]
        #            if self.u.dof[J] > 0.0:
        #                print 'calc qu j= ',j,' u.dof[',J,']= ',self.u.dof[J]
        #            manu[eN,k]+= self.u.dof[J]*self.q['v'][eN,k,j]
        #        #end j
        #        print """calc q[u][%d,%d]= %g """ % (eN,k,manu[eN,k])
        #    #end k
        # end eN
        # mwf debug end
    def calculateFiniteElementU_elementBoundaries(self):
        """
        Calculate u and grad(u)
        """
        self.u.getValuesTrace(self.ebq['v'],self.ebq['u'])
        self.u.getGradientValuesTrace(self.ebq['grad(v)'],self.ebq['grad(u)'])
    def calculateFiniteElementPhi(self):
        """
        Calculate phi, grad(phi) and the (outer) product grad(phi) x grad(w)
        """
        self.phi.getValues(self.q['v'],
                           self.q['phi'])
        self.phi.getGradientValues(self.q['grad(v)'],
                                   self.q['grad(phi)'])
        self.phi.getGradientTensorValues(self.q['grad(v)*grad(w)'],
                                         self.q['grad(phi)*grad(w)'])
    def calculateFiniteElementPhi_elementBoundaries(self):
        """
        Calculate phi, grad(phi) and the (outer) product grad(phi) x grad(w)
        """
        self.phi.getValuesTrace(self.ebq['v'],self.ebq['phi'])
        self.phi.getGradientValuesTrace(self.ebq['grad(v)'],self.ebq['grad(phi)'])
    def updateMass(self,residual):
        """
        Calculate the element mass accumulation integral and update
        the element residual
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateMass(#self.mesh.nElements_global,
                                #self.nQuadraturePoints_element,
                                #self.nDOF_element,
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
                #mwf debug
                #mwf added
                #print 'q[mt]= ',self.q['mt']
                #print 'q[dx_m]= ',self.q['dx_m']
                #print 'q[mt*dx_m]= ',self.q['mt*dx_m']
                #print 'q[w]= ',self.q['w']
                #print 'restmp= ',self.restemp
                #print 'residual= ',residual
                print 'max(restemp-residual)= ',max(self.restemp.flat[:]-residual.flat[:])
                if max(self.restemp.flat[:]-residual.flat[:]) > 1.0e-6:
                       print 'q[mt]= ',self.q['mt']
                       print 'q[w]= ',self.q['w']
                       print 'q[mt*dx_m]= ',self.q['mt*dx_m']
                       print 'DT = ',self.timeIntegration.DT
                print 'restemp-residual= '
                #mwf end added
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
    def updateAdvection(self,residual):
        """
        Calculate the element advection integral and update element
        residual
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateAdvection(self.mesh.nElements_global,
                                     self.nQuadraturePoints_element,
                                     self.nDOF_element,
                                     self.nSpace_global,
                                     self.q['f*dx_f'],
                                     self.q['grad(w)'],
                                     residual)
        #mwf debug
        #print 'updateAdvection grad(w)[eN,k,i,d] '
        #for eN in range(self.mesh.nElements_global):
        #    for k in range(self.nQuadraturePoints_element):
        #        for i in range(self.nDOF_element):
        #                for d in range(self.nSpace_global):
        #                    print """qrad(w)[%d,%d,%d,%d]= %g
        #                    """ % (eN,k,i,d,self.q['grad(w)'][eN,k,i,d])
        #end eN
        if self.testUpdateIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        self.restemp[eN,i] -= (
                            Numeric.dot(self.q['f*dx_f'][eN,k],
                                        self.q['grad(w)'][eN,k,i]))
            if min(self.restemp.flat == residual.flat) == 0:
                #mwf added
                print "advection error"
                print 'max(advrestemp-residual)= ',max(self.restemp.flat[:]-residual.flat[:])
                print 'restemp-residual= '
                print self.restemp - residual
                print "advection error"
    def updateAdvectionJacobian(self,jacobian):
        """
        Calculate the jacobian of the element advection integral and
        update element jacobian
        """
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp.flat[:]=jacobian.flat
        femIntegrals.updateAdvectionJacobian(self.mesh.nElements_global,
                                             self.nQuadraturePoints_element,
                                             self.nDOF_element,
                                             self.nSpace_global,
                                             self.q['df*dx_f'],
                                             self.q['v*grad(w)'],
                                             jacobian)
        if self.testUpdateJacobianIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for j in range(self.nDOF_element):
                            self.jacobiantemp[eN,i,j] -= (
                                Numeric.dot(self.q['df*dx_f'][eN,k],
                                            self.q['v*grad(w)'][eN,k,j,i]))
            if min(self.jacobiantemp.flat == jacobian.flat) == 0:
                print "advection jacobian error"
                print self.jacobiantemp - jacobian
                print "advection jacobian error"
    def updateDiffusion(self,residual):
        """
        Calculate the element diffusion integral and update the
        element residual
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateDiffusion(self.mesh.nElements_global,
                                     self.nQuadraturePoints_element,
                                     self.nDOF_element,
                                     self.nSpace_global,
                                     self.q['a*dx_a'],
                                     self.q['grad(phi)*grad(w)'],
                                     residual)
        if self.testUpdateIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
#                          self.restemp[eN,i] += (
#                              Numeric.innerproduct(self.q['a*dx_a'][eN,k].flat,
#                                                   self.q['grad(phi)*grad(w)'][eN,k,i].flat))
                        for I in range(self.nSpace_global):
                            for J in range(self.nSpace_global):
                                self.restemp[eN,i] += (
                                    self.q['a*dx_a'][eN,k,I,J]*
                                    self.q['grad(phi)*grad(w)'][eN,k,i,J,I])
            if min(self.restemp.flat == residual.flat) == 0:
                print "diffusion error"
                print self.restemp - residual
                print "diffusion error"
    def updateDiffusionJacobian(self,jacobian):
        """
        Calculate the jacobian of the element diffusion integral
        and update element jacobian
        """
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp.flat[:]=jacobian.flat
        femIntegrals.updateDiffusionJacobian(self.mesh.nElements_global,
                                             self.nQuadraturePoints_element,
                                             self.nDOF_element,
                                             self.nSpace_global,
                                             self.phi.femSpace.dofMap.l2g,
                                             self.q['a*dx_a'],
                                             self.q['da*dx_a'],
                                             self.q['grad(phi)*grad(w)'],
                                             self.dphi.dof,
                                             self.q['v'],
                                             self.q['grad(v)*grad(w)'],
                                             jacobian)
        if self.testUpdateJacobianIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for j in range(self.nDOF_element_phi):
                            J = self.trialSpace.dofMap.l2g[eN,j]
                            #print "dphiProductPython"
                            #print Numeric.innerproduct(
                            #    self.q['a*dx_a'][eN,k].flat,
                            #    self.q['grad(v)*grad(w)'][eN,k,j,i].flat)
                            self.jacobiantemp[eN,i,j] += (
                                Numeric.innerproduct(
                                self.q['da*dx_a'][eN,k].flat,
                                self.q['grad(phi)*grad(w)'][eN,k,i].flat)*self.q['v'][eN,k,j] + 
                                self.dphi.dof[J]*Numeric.innerproduct(
                                self.q['a*dx_a'][eN,k].flat,
                                self.q['grad(v)*grad(w)'][eN,k,j,i].flat))
            if min(self.jacobiantemp.flat == jacobian.flat) == 0:
                print "diffusion jacobian error"
                print self.jacobiantemp - jacobian
                print "diffusion jacobian error"
    def updateReaction(self,residual):
        """
        Calculate the element reaction integral and the update element
        residual
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateReaction(self.mesh.nElements_global,
                                    self.nQuadraturePoints_element,
                                    self.nDOF_element,
                                    self.q['r*dx_r'],
                                    self.q['w'],
                                    residual)
        if self.testUpdateIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        self.restemp[eN,i] += (self.q['r*dx_r'][eN,k]*
                                           self.q['w'][eN,k,i])
            if min(self.restemp.flat == residual.flat) == 0:
                print "reaction error"
                print self.restemp - residual
                print "reaction error"
    def updateReactionJacobian(self,jacobian):
        """
        Calculate the jacobian of the element reaction integral and
        the update element jacobian
        """
        if self.testUpdateJacobianIntegrals:
            self.jacobiantemp.flat[:]=jacobian.flat
        femIntegrals.updateReactionJacobian(self.mesh.nElements_global,
                                            self.nQuadraturePoints_element,
                                            self.nDOF_element,
                                            self.q['dr*dx_r'],
                                            self.q['v*w'],
                                            jacobian)
        if self.testUpdateJacobianIntegrals:
            for eN in range(self.mesh.nElements_global):
                for i in range(self.nDOF_element):
                    for k in range(self.nQuadraturePoints_element):
                        for j in range(self.nDOF_element):
                            self.jacobiantemp[eN,i,j] += (self.q['dr*dx_r'][eN,k]*
                                                          self.q['v*w'][eN,k,j,i])
            if min(self.jacobiantemp.flat == jacobian.flat) == 0:
                print "reaction jacobian error"
                print self.jacobiantemp - jacobian
                print "reaction jacobian error"
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
            self.restemp.flat[:]=residual.flat
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
    def updateAdvectionInteriorElementBoundaryFlux(self,residual):
        """
        Update the residual with the flux through the element boundaries
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateInteriorElementBoundaryFlux(self.mesh.nInteriorElementBoundaries_global,
                                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                       self.mesh.nElementBoundaries_element,
                                                       self.nDOF_element,
                                                       self.mesh.interiorElementBoundariesArray,
                                                       self.mesh.elementBoundaryElementsArray,
                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                       self.ebq_global['advectiveFlux*dx_f'],
                                                       self.ebq['w'],
                                                       residual)
        if self.testUpdateIntegrals:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                for i in range(self.nDOF_element):
                    for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                        self.restemp[left_eN_global,i] += self.ebq_global['advectiveFlux*dx_f'][ebN,k]*self.ebq['w'][left_eN_global,left_ebN_element,k,i]
                        self.restemp[right_eN_global,i] -= self.ebq_global['advectiveFlux*dx_f'][ebN,k]*self.ebq['w'][right_eN_global,right_ebN_element,k,i]
            if min(self.restemp.flat == residual.flat) == 0:
                print "advection interior element boundary flux error"
                print self.restemp - residual
                print "advection interior element boundary flux error"
    def updateAdvectionExteriorElementBoundaryFlux(self,residual):
        """
        Update the residual with the flux through the element boundaries
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateExteriorElementBoundaryFlux(self.mesh.nExteriorElementBoundaries_global,
                                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                       self.mesh.nElementBoundaries_element,
                                                       self.nDOF_element,
                                                       self.mesh.exteriorElementBoundariesArray,
                                                       self.mesh.elementBoundaryElementsArray,
                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                       self.ebq_global['advectiveFlux*dx_f'],
                                                       self.ebq['w'],
                                                       residual)
        if self.testUpdateIntegrals:
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                for i in range(self.nDOF_element):
                    for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                        self.restemp[eN_global,i] += self.ebq_global['advectiveFlux*dx_f'][ebN,k]*self.ebq['w'][eN_global,ebN_element,k,i]
            if min(self.restemp.flat == residual.flat) == 0:
                print "advection exterior element boundary flux error"
                print self.restemp - residual
                print "advection exterior element boundary flux error"
    def updateDiffusionInteriorElementBoundaryFlux(self,residual):
        """
        Update the residual with the flux through the element boundaries
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateInteriorElementBoundaryFlux(self.mesh.nInteriorElementBoundaries_global,
                                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                       self.mesh.nElementBoundaries_element,
                                                       self.nDOF_element,
                                                       self.mesh.interiorElementBoundariesArray,
                                                       self.mesh.elementBoundaryElementsArray,
                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                       self.ebq_global['diffusiveFlux*dx_a'],
                                                       self.ebq['w'],
                                                       residual)
        if self.testUpdateIntegrals:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                for i in range(self.nDOF_element):
                    for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                        self.restemp[left_eN_global,i]  += self.ebq_global['diffusiveFlux*dx_a'][ebN,k]*self.ebq['w'][left_eN_global,left_ebN_element,k,i]
                        self.restemp[right_eN_global,i] -= self.ebq_global['diffusiveFlux*dx_a'][ebN,k]*self.ebq['w'][right_eN_global,right_ebN_element,k,i]
            if min(self.restemp.flat == residual.flat) == 0:
                print "diffusion interior element boundary flux error"
                print self.restemp - residual
                print "diffusion interior element boundary flux error"
    def updateDiffusionExteriorElementBoundaryFlux(self,residual):
        """
        Update the residual with the flux through the element boundaries
        """
        if self.testUpdateIntegrals:
            self.restemp.flat[:]=residual.flat
        femIntegrals.updateExteriorElementBoundaryFlux(self.mesh.nExteriorElementBoundaries_global,
                                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                       self.mesh.nElementBoundaries_element,
                                                       self.nDOF_element,
                                                       self.mesh.exteriorElementBoundariesArray,
                                                       self.mesh.elementBoundaryElementsArray,
                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                       self.ebq_global['diffusiveFlux*dx_a'],
                                                       self.ebq['w'],
                                                       residual)
        if self.testUpdateIntegrals:
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                for i in range(self.nDOF_element):
                    for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                        self.restemp[eN_global,i] += self.ebq_global['diffusiveFlux*dx_a'][ebN,k]*self.ebq['w'][eN_global,ebN_element,k,i]
            if min(self.restemp.flat == residual.flat) == 0:
                print "diffusion exterior element boundary flux error"
                print self.restemp - residual
                print "diffusion exterior element boundary flux error"
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
    def getFlowVelocity(self):
#          for eN in range(self.mesh.nElements_global):
#              for  k  in range(self.nQuadraturePoints_element):
#                  for  I in range(self.nSpace_global):
#                      self.q['velocity'][eN,k,I] = self.q['f'][eN,k,I]
#                      for  J in range(self.nSpace_global):
#                          self.q['velocity'][eN,k,I] -= (
#                              self.q['a'][eN,k,I,J]*self.q['grad(phi)'][eN,k,J])
        femIntegrals.calculateFlowVelocity(self.mesh.nElements_global,
                                           self.nQuadraturePoints_element,
                                           self.nSpace_global,
                                           self.q['f'],
                                           self.q['a'],
                                           self.q['grad(phi)'],
                                           self.q['velocity'])
        return self.q['velocity']
    def getFlowVelocityElementBoundaries(self):
        self.updateElementBoundaryCoefficients()
        for eN in range(self.mesh.nElements_global):
            for ebN in range(self.mesh.nElementBoundaries_element):
                for  k  in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    for  I in range(self.nSpace_global):
                        self.ebq['velocity'][eN,ebN,k,I] = self.ebq['f'][eN,ebN,k,I]
                        for  J in range(self.nSpace_global):
                            self.ebq['velocity'][eN,ebN,k,I] -= (
                                self.ebq['a'][eN,ebN,k,I,J]*self.ebq['grad(phi)'][eN,ebN,k,J])
        return self.ebq['velocity']
    def calculateInteriorNumericalAdvectiveFlux(self):
        femIntegrals.calculateInteriorNumericalAdvectiveFlux(self.mesh.nInteriorElementBoundaries_global,
                                                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                             self.mesh.nElementBoundaries_element,
                                                             self.nSpace_global,
                                                             self.mesh.interiorElementBoundariesArray,
                                                             self.mesh.elementBoundaryElementsArray,
                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                             self.ebq['n'],
                                                             self.ebq['f'],
                                                             self.ebq['df'],
                                                             self.ebq_global['dx_f'],
                                                             self.ebq_global['advectiveFlux*dx_f'])
        if self.testNumericalFlux:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                left_normals   = self.ebq['n'][left_eN_global,left_ebN_element]
                right_normals  = self.ebq['n'][right_eN_global,right_ebN_element]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    left_char_speed = dot(self.ebq['df'][left_eN_global,left_ebN_element,k],left_normals[k])
                    right_char_speed = dot(self.ebq['df'][right_eN_global,right_ebN_element,k],left_normals[k])
                    speed=right_char_speed
                    if max(abs(left_char_speed),abs(right_char_speed)) == abs(left_char_speed):
                        speed=left_char_speed
                    if speed >= 0.0:
                        self.fluxtemp[ebN,k] = (dot(self.ebq['f'][left_eN_global,left_ebN_element,k],
                                                    left_normals[k])
                                                *
                                                self.ebq_global['dx_f'][ebN,k])
                    else:
                        self.fluxtemp[ebN,k] = (dot(self.ebq['f'][right_eN_global,right_ebN_element,k],
                                                    left_normals[k])
                                                *
                                                self.ebq_global['dx_f'][ebN,k])
                    if self.fluxtemp[ebN,k] != self.ebq_global['advectiveFlux*dx_f'][ebN,k]:
                        print "interior numerical advective flux error"
                        print self.fluxtemp[ebN,k] - self.ebq_global['advectiveFlux*dx_f'][ebN,k]
                        print "interior numerical advective flux error"
    def updateInteriorNumericalAdvectiveFluxJacobian(self,fluxJacobian):
        if self.testNumericalFluxJacobian:
            self.fluxJacobiantemp.flat[:]=fluxJacobian.flat
        femIntegrals.updateInteriorNumericalAdvectiveFluxJacobian(self.mesh.nInteriorElementBoundaries_global,
                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                  self.mesh.nElementBoundaries_element,
                                                                  self.nDOF_element,
                                                                  self.nSpace_global,
                                                                  self.mesh.interiorElementBoundariesArray,
                                                                  self.mesh.elementBoundaryElementsArray,
                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                  self.ebq['n'],
                                                                  self.ebq['df'],
								  self.ebq['v'],
                                                                  self.ebq_global['dx_f'],
                                                                  fluxJacobian)
        if self.testNumericalFluxJacobian:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                left_normals   = self.ebq['n'][left_eN_global,left_ebN_element]
                right_normals  = self.ebq['n'][right_eN_global,right_ebN_element]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    left_char_speed = dot(self.ebq['df'][left_eN_global,left_ebN_element,k],left_normals[k])
                    right_char_speed = dot(self.ebq['df'][right_eN_global,right_ebN_element,k],left_normals[k])
                    speed=right_char_speed
                    if max(math.fabs(left_char_speed),math.fabs(right_char_speed)) == math.fabs(left_char_speed):
                        speed=left_char_speed
                    if speed >= 0.0:
                        for j in range(self.nDOF_element):
                            self.fluxJacobiantemp[ebN,0,k,j] += (dot(self.ebq['df'][left_eN_global,left_ebN_element,k],
                                                                    left_normals[k])
                                                                *
                                                                self.ebq_global['dx_f'][ebN,k]
                                                                *
                                                                self.ebq['v'][left_eN_global,left_ebN_element,k,j])
                    else:
                        for j in range(self.nDOF_element):
                             self.fluxJacobiantemp[ebN,1,k,j] += (dot(self.ebq['df'][right_eN_global,right_ebN_element,k],
                                                                     left_normals[k])
                                                                 *
                                                                 self.ebq_global['dx_f'][ebN,k]
                                                                 *
                                                                 self.ebq['v'][right_eN_global,right_ebN_element,k,j])
                    if fluxJacobian[ebN,0,k,j] != self.fluxJacobiantemp[ebN,0,k,j]:
                        print "left interior advective flux jacobian error"
                        print fluxJacobian[ebN,0,k,j] - self.fluxJacobiantemp[ebN,0,k,j]
                        print "left interior advective flux jacobian error"
                    if fluxJacobian[ebN,1,k,j] != self.fluxJacobiantemp[ebN,1,k,j]:
                        print "right interior advective flux jacobian error"
                        print fluxJacobian[ebN,1,k,j] - self.fluxJacobiantemp[ebN,1,k,j]
                        print "right interior advective flux jacobian error"
    def calculateExteriorNumericalAdvectiveFlux(self):
#          if self.fluxBoundaryConditions == 'noFlow':
#              for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#                  ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#                  eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  normals   = self.ebq['n'][eN_global,ebN_element]
#                  for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                      self.ebq_global['advectiveFlux*dx_f'][ebN,k] = 0.0
        femIntegrals.calculateExteriorNumericalAdvectiveFlux(self.mesh.nExteriorElementBoundaries_global,
                                                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                             self.mesh.nElementBoundaries_element,
                                                             self.nSpace_global,
                                                             self.mesh.exteriorElementBoundariesArray,
                                                             self.mesh.elementBoundaryElementsArray,
                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                             self.inflowBoundary,
                                                             self.inflowFlux,
                                                             self.ebq['n'],
                                                             self.ebq['f'],
                                                             self.ebq['df'],
                                                             self.ebq_global['dx_f'],
                                                             self.ebq_global['advectiveFlux*dx_f'])
    def setInflowFlux(self):
        femIntegrals.setInflowFlux(self.mesh.nExteriorElementBoundaries_global,
                                   self.nElementBoundaryQuadraturePoints_elementBoundary,
                                   self.mesh.exteriorElementBoundariesArray,
                                   self.ebq_global['advectiveFlux*dx_f'],
                                   self.inflowFlux)
    def updateExteriorNumericalAdvectiveFluxJacobian(self,fluxJacobian):
#          if self.fluxBoundaryConditions == 'noFlow':
#              for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
#                  ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
#                  eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
#                  ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
#                  normals   = self.ebq['n'][eN_global,ebN_element]
#                  for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
#                      for j in range(self.nDOF_element):
#                          fluxJacobian[ebN,0,k,j] = 0.0
        if self.fluxBoundaryConditions == 'outFlow':
            if self.testNumericalFluxJacobian:
                self.fluxJacobiantemp.flat[:]=fluxJacobian.flat
            femIntegrals.updateExteriorNumericalAdvectiveFluxJacobian(self.mesh.nExteriorElementBoundaries_global,
                                                                      self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                      self.mesh.nElementBoundaries_element,
                                                                      self.nDOF_element,
                                                                      self.nSpace_global,
                                                                      self.mesh.exteriorElementBoundariesArray,
                                                                      self.mesh.elementBoundaryElementsArray,
                                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                      self.inflowBoundary,
                                                                      self.ebq['n'],
                                                                      self.ebq['df'],
								      self.ebq['v'],
                                                                      self.ebq_global['dx_f'],
                                                                      fluxJacobian)
            if  self.testNumericalFluxJacobian:
                for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                    ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                    eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                    ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                    normals   = self.ebq['n'][eN_global,ebN_element]
                    for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                        for j in range(self.nDOF_element):
                            self.fluxJacobiantemp[ebN,0,k,j] += (dot(self.ebq['df'][eN_global,ebN_element,k],
                                                                     normals[k])
                                                                 *
                                                                 self.ebq_global['dx_f'][ebN,k]
                                                                 *
                                                                 self.ebq['v'][eN_global,ebN_element,k,j])
                    if self.fluxJacobiantemp[ebN,0,k,j] != fluxJacobian[ebN,0,k,j]:
                        print  "error in exterior advective flux jacobian"
                        print self.fluxJacobiantemp[ebN,0,k,j] - fluxJacobian[ebN,0,k,j]
                        print  "error in exterior advective flux jacobian"
    def calculateInteriorNumericalDiffusiveFlux(self):
        femIntegrals.calculateInteriorNumericalDiffusiveFlux(self.mesh.nInteriorElementBoundaries_global,
                                                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                             self.mesh.nElementBoundaries_element,
                                                             self.nSpace_global,
                                                             self.mesh.interiorElementBoundariesArray,
                                                             self.mesh.elementBoundaryElementsArray,
                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                             self.ebq['n'],
                                                             self.ebq['a'],
                                                             self.ebq['grad(phi)'],
                                                             self.ebq['u'],
                                                             self.ebq_global['penalty'],
                                                             self.ebq_global['dx_a'],
                                                             self.ebq_global['diffusiveFlux*dx_a'])
        if self.testNumericalFlux:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                left_normals   = self.ebq['n'][left_eN_global,left_ebN_element]
                right_normals  = self.ebq['n'][right_eN_global,right_ebN_element]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.fluxtemp[ebN,k] = (-0.5*dot(tenvec(self.ebq['a'][left_eN_global,left_ebN_element,k],
                                                            self.ebq['grad(phi)'][left_eN_global,left_ebN_element,k])
                                                     +
                                                     tenvec(self.ebq['a'][right_eN_global,right_ebN_element,k],
                                                            self.ebq['grad(phi)'][right_eN_global,right_ebN_element,k]),
                                                     left_normals[k])*self.ebq_global['dx_a'][ebN,k]
                                            +(2.0/self.eb_global['heb'][ebN])*(self.ebq['u'][left_eN_global,left_ebN_element,k]-
                                                                               self.ebq['u'][right_eN_global,right_ebN_element,k])*self.ebq_global['dx_a'][ebN,k])
                if self.fluxtemp[ebN,k] != self.ebq_global['diffusiveFlux*dx_a'][ebN,k]:
                    print "error  in interior numerical diffusive  flux"
                    print self.fluxtemp[ebN,k] - self.ebq_global['diffusiveFlux*dx_a'][ebN,k]
                    print "error  in interior numerical diffusive  flux"
    def updateInteriorNumericalDiffusiveFluxJacobian(self,fluxJacobian):
        if self.testNumericalFluxJacobian:
            self.fluxJacobiantemp.flat[:]=fluxJacobian.flat
        femIntegrals.updateInteriorNumericalDiffusiveFluxJacobian(self.mesh.nInteriorElementBoundaries_global,
                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                  self.mesh.nElementBoundaries_element,
                                                                  self.nDOF_element,
                                                                  self.nSpace_global,
                                                                  self.u.femSpace.dofMap.l2g,
                                                                  self.mesh.interiorElementBoundariesArray,
                                                                  self.mesh.elementBoundaryElementsArray,
                                                                  self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                  self.ebq['n'],
                                                                  self.ebq['a'],
                                                                  self.ebq['da'],
                                                                  self.ebq['grad(phi)'],
                                                                  self.dphi.dof,
                                                                  self.ebq['v'],
                                                                  self.ebq['grad(v)'],
                                                                  self.ebq_global['penalty'],
                                                                  self.ebq_global['dx_a'],
                                                                  fluxJacobian)
        if self.testNumericalFluxJacobian:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                left_normals   = self.ebq['n'][left_eN_global,left_ebN_element]
                right_normals  = self.ebq['n'][right_eN_global,right_ebN_element]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    for j in range(self.nDOF_element):
                        left_J = self.trialSpace.dofMap.l2g[left_eN_global,j]
                        right_J = self.trialSpace.dofMap.l2g[right_eN_global,j]
 #                         self.fluxJacobiantemp[ebN,0,k,j] += (-0.5*dot(tenvec(self.ebq['da'][left_eN_global,left_ebN_element,k],
#                                                                               self.ebq['grad(phi)'][left_eN_global,left_ebN_element,k])*self.ebq['v'][left_eN_global,left_ebN_element,k,j]+
#                                                                        self.dphi.dof[left_J]*tenvec(self.ebq['a'][left_eN_global,left_ebN_element,k],
#                                                                                                     self.ebq['grad(v)'][left_eN_global,left_ebN_element,k,j]),
#                                                                        left_normals[k])*self.ebq_global['dx_a'][ebN,k]
#                                                               +(2.0/self.eb_global['heb'][ebN])*self.ebq['v'][left_eN_global,left_ebN_element,k,j]*self.ebq_global['dx_a'][ebN,k])
#                          self.fluxJacobiantemp[ebN,1,k,j] += (-0.5*dot(tenvec(self.ebq['da'][right_eN_global,right_ebN_element,k],
#                                                                               self.ebq['grad(phi)'][right_eN_global,right_ebN_element,k])*self.ebq['v'][right_eN_global,right_ebN_element,k,j]+
#                                                                        self.dphi.dof[right_J]*tenvec(self.ebq['a'][right_eN_global,right_ebN_element,k],
#                                                                                                      self.ebq['grad(v)'][right_eN_global,right_ebN_element,k,j]),
#                                                                        left_normals[k])*self.ebq_global['dx_a'][ebN,k]
#                                                               -(2.0/self.eb_global['heb'][ebN])*self.ebq['v'][right_eN_global,right_ebN_element,k,j]*self.ebq_global['dx_a'][ebN,k])
                        self.fluxJacobiantemp[ebN,0,k,j] += (-0.5*dot(tenvec(self.ebq['da'][left_eN_global,left_ebN_element,k],
                                                                             self.ebq['grad(phi)'][left_eN_global,left_ebN_element,k])*self.ebq['v'][left_eN_global,left_ebN_element,k,j]+
                                                                      self.dphi.dof[left_J]*tenvec(self.ebq['a'][left_eN_global,left_ebN_element,k],
                                                                                                   self.ebq['grad(v)'][left_eN_global,left_ebN_element,k,j]),
                                                                      left_normals[k])*self.ebq_global['dx_a'][ebN,k]
                                                             +self.ebq_global['penalty'][ebN,k]*self.ebq['v'][left_eN_global,left_ebN_element,k,j]*self.ebq_global['dx_a'][ebN,k])
                        self.fluxJacobiantemp[ebN,1,k,j] += (-0.5*dot(tenvec(self.ebq['da'][right_eN_global,right_ebN_element,k],
                                                                             self.ebq['grad(phi)'][right_eN_global,right_ebN_element,k])*self.ebq['v'][right_eN_global,right_ebN_element,k,j]+
                                                                      self.dphi.dof[right_J]*tenvec(self.ebq['a'][right_eN_global,right_ebN_element,k],
                                                                                                    self.ebq['grad(v)'][right_eN_global,right_ebN_element,k,j]),
                                                                      left_normals[k])*self.ebq_global['dx_a'][ebN,k]
                                                             -self.ebq_global['penalty'][ebN,k]*self.ebq['v'][right_eN_global,right_ebN_element,k,j]*self.ebq_global['dx_a'][ebN,k])
                    if self.fluxJacobiantemp[ebN,0,k,j] != fluxJacobian[ebN,0,k,j]:
                        print "error in left interior numerical diffusive flux jacobian",ebN,0,k,j
                        print self.fluxJacobiantemp[ebN,0,k,j] - fluxJacobian[ebN,0,k,j]
                        print "error in left interior numerical diffusive flux jacobian"
                    if self.fluxJacobiantemp[ebN,1,k,j] != fluxJacobian[ebN,1,k,j]:
                        print "error in right interior numerical diffusive flux jacobian",ebN,1,k,j
                        print self.fluxJacobiantemp[ebN,1,k,j] - fluxJacobian[ebN,1,k,j]
                        print "error in right interior numerical diffusive flux jacobian"
    def calculateExteriorNumericalDiffusiveFlux(self):
        pass
        #TODO
#          femIntegrals.calculateExteriorNumericalDiffusiveFlux(self.mesh.nExteriorElementBoundaries_global,
#                                                               self.nElementBoundaryQuadraturePoints_elementBoundary,
#                                                               self.mesh.nElementBoundaries_element,
#                                                               self.nSpace_global,
#                                                               self.mesh.exteriorElementBoundariesArray,
#                                                               self.mesh.elementBoundaryElementsArray,
#                                                               self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                               self.ebq['n'],
#                                                               self.ebq['a'],
#                                                               self.ebq['grad(phi)'],
#                                                               self.ebq['u'],
#                                                               self.ebq_global['penalty'],
#                                                               self.ebq_global['dx_a'], 
#                                                               self.ebq_global['diffusiveFlux*dx_a'])
        #if self.fluxBoundaryConditions == 'noFlow':
        #    for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
        #        ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
        #        eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
        #        ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
        #        normals   = self.ebq['n'][eN_global,ebN_element]
        #        for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
        #            self.ebq_global['diffusiveFlux*dx_a'][ebN,k] = 0.0
        #if self.fluxBoundaryConditions == 'outFlow':
        #    for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
        #        ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
        #        eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
        #        ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
        #        normals   = self.ebq['n'][eN_global,ebN_element]
        #        for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
        #            self.ebq_global['diffusiveFlux*dx_a'][ebN,k] = 0.0
    def updateExteriorNumericalDiffusiveFluxJacobian(self,fluxJacobian):
        pass
        #TODO
#          femIntegrals.updateExteriorNumericalDiffusiveFluxJacobian(self.mesh.nExteriorElementBoundaries_global,
#                                                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
#                                                                       self.mesh.nElementBoundaries_element,
#                                                                       self.nDOF_element,
#                                                                       self.nSpace_global,
#                                                                       self.u.femSpace.dofMap.l2g,
#                                                                       self.mesh.exteriorElementBoundariesArray,
#                                                                       self.mesh.elementBoundaryElementsArray,
#                                                                       self.mesh.elementBoundaryLocalElementBoundariesArray,
#                                                                       self.ebq['n'],
#                                                                       self.ebq['a'],
#                                                                       self.ebq['da'],
#                                                                       self.ebq['grad(phi)'],
#                                                                       self.dphi.dof,
#                                                                       self.ebq['v'],
#                                                                       self.ebq['grad(v)'],
#                                                                       self.ebq_global['penalty'],
#                                                                       self.ebq_global['dx_a'],
#                                                                       fluxJacobian)
        #if self.fluxBoundaryConditions == 'noFlow':
        #    for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
        #        ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
        #        eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
        #        ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
        #        normals   = self.ebq['n'][eN_global,ebN_element]
        #        for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
        #            for j in range(self.nDOF_element):
        #                fluxJacobian[ebN,0,k,j] += 0.0
        #if self.fluxBoundaryConditions == 'outFlow':
        #    for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
        #        ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
        #        eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
        #        ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
        #        normals   = self.ebq['n'][eN_global,ebN_element]
        #        for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
        #            for j in range(self.nDOF_element):
        #                fluxJacobian[ebN,0,k,j] += 0.0
    def getConservationResidualPWC(self):
        self.updateElementBoundaryCoefficients()
        self.ebq_global['conservationVelocityPWC'].flat[:]=0.0
        self.eb_global['conservationFluxPWC'].flat[:]=0.0
	self.ebq['velocity'].flat[:] = 0.0
        for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
            ebN = self.mesh.interiorElementBoundariesArray[ebNI]
            left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
            left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
            left_normals   = self.ebq['n'][left_eN_global,left_ebN_element]
            right_normals  = self.ebq['n'][right_eN_global,right_ebN_element]
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                self.ebq_global['conservationVelocityPWC'][ebN,k,:] = (self.ebq_global['velocityAverage'][ebN,k] -
                                                                       (self.e['conservationCorrectionPWC'][left_eN_global]-
                                                                        self.e['conservationCorrectionPWC'][right_eN_global])*
                                                                       left_normals[k])
		self.ebq['velocity'][left_eN_global,left_ebN_element,k,:] = self.ebq_global['conservationVelocityPWC'][ebN,k,:]
		self.ebq['velocity'][right_eN_global,right_ebN_element,k,:] = self.ebq_global['conservationVelocityPWC'][ebN,k,:]	    
                self.eb_global['conservationFluxPWC'][ebN] += dot(self.ebq_global['conservationVelocityPWC'][ebN,k],left_normals[k])*self.ebq_global['dx_f'][ebN,k]
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            normals   = self.ebq['n'][eN_global,ebN_element]
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                self.ebq_global['conservationVelocityPWC'][ebN,k,:] = (self.ebq_global['velocityAverage'][ebN,k] -
                                                                       self.e['conservationCorrectionPWC'][eN_global]*normals[k])
		self.ebq['velocity'][eN_global,ebN_element,k,:] = self.ebq_global['conservationVelocityPWC'][ebN,k,:]
                self.eb_global['conservationFluxPWC'][ebN] += dot(self.ebq_global['conservationVelocityPWC'][ebN,k],normals[k])*self.ebq_global['dx_f'][ebN,k]
        self.e['conservationResidual'].flat[:]=0.0
        for eN in range(self.mesh.nElements_global):
            for i in range(self.nDOF_element):
                self.e['conservationResidual'][eN]+=self.elementResidual[eN,i]
        for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
            ebN = self.mesh.interiorElementBoundariesArray[ebNI]
            left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
            left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
            self.e['conservationResidual'][left_eN_global]+=self.eb_global['conservationFluxPWC'][ebN]
            self.e['conservationResidual'][right_eN_global]-=self.eb_global['conservationFluxPWC'][ebN]
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            self.e['conservationResidual'][eN_global]+=self.eb_global['conservationFluxPWC'][ebN]
    def getConservationJacobianPWC(self):
        #initialize matrix
        for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
            ebN = self.mesh.interiorElementBoundariesArray[ebNI]
            left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
            left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
            self.conservationJacobianPWC[left_eN_global,left_eN_global] = 0.0
            self.conservationJacobianPWC[left_eN_global,right_eN_global] = 0.0
            self.conservationJacobianPWC[right_eN_global,left_eN_global] = 0.0
            self.conservationJacobianPWC[right_eN_global,right_eN_global] = 0.0
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            self.conservationJacobianPWC[eN_global,eN_global] = 0.0
        #set matrix
        for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
            ebN = self.mesh.interiorElementBoundariesArray[ebNI]
            left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
            left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
            w=0.0
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                w+=self.ebq_global['dx_f'][ebN,k]
            self.conservationJacobianPWC[left_eN_global,left_eN_global] -= w#1.0
            self.conservationJacobianPWC[left_eN_global,right_eN_global] += w#1.0
            self.conservationJacobianPWC[right_eN_global,left_eN_global] += w#1.0
            self.conservationJacobianPWC[right_eN_global,right_eN_global] -= w#1.0
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            w=0.0
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                w+=self.ebq_global['dx_f'][ebN,k]
            self.conservationJacobianPWC[eN_global,eN_global] -= w#1.0
        self.conservationSolver.prepare()
    def getConservationResidualPWL(self):
	femIntegrals.calculateConservationResidualPWL(self.mesh.nElements_global,
						      self.mesh.nInteriorElementBoundaries_global,
						      self.mesh.nExteriorElementBoundaries_global,
						      self.nElementBoundaryQuadraturePoints_elementBoundary,
						      self.mesh.nElementBoundaries_element,
						      self.mesh.nNodes_element,
						      self.nSpace_global,
						      self.mesh.interiorElementBoundariesArray,
						      self.mesh.exteriorElementBoundariesArray,
						      self.mesh.elementBoundaryElementsArray,
						      self.mesh.elementBoundaryLocalElementBoundariesArray,
						      self.mesh.elementNodesArray,
						      self.nodeStarElementsArray,
						      self.nodeStarElementNeighborsArray,
                                                      self.nodeStarOffset,
						      self.mesh.nElements_node,
						      self.elementResidual,
						      self.ebq_global['velocityAverage'],
						      self.starU,
						      self.ebq['w'],
						      self.ebq_global['n'],
						      self.ebq_global['dx_f'],
						      self.e['conservationResidual'],
						      self.starR,
						      self.ebq_global['conservationVelocityPWL'],
						      self.ebq['velocity'])
        self.ebq_global['velocity'] = self.ebq_global['conservationVelocityPWL']
    def getConservationJacobianPWL(self):
	femIntegrals.calculateConservationJacobianPWL(self.mesh.nNodes_global,
                                                      self.nNodes_internal,
                                                      self.mesh.nElements_global,
						      self.mesh.nInteriorElementBoundaries_global,
						      self.mesh.nExteriorElementBoundaries_global,
						      self.nElementBoundaryQuadraturePoints_elementBoundary,
						      self.mesh.nElementBoundaries_element,
						      self.mesh.nNodes_element,
						      self.nSpace_global,
						      self.mesh.interiorElementBoundariesArray,
						      self.mesh.exteriorElementBoundariesArray,
						      self.mesh.elementBoundaryElementsArray,
						      self.mesh.elementBoundaryLocalElementBoundariesArray,
						      self.mesh.elementNodesArray,
						      self.nodeStarElementsArray,
						      self.nodeStarElementNeighborsArray,
                                                      self.nodeStarOffset,
                                                      self.nodeStarJacobianOffset,
						      self.mesh.nElements_node,
                                                      self.internalNodesArray,
						      self.ebq['w'],
						      self.ebq_global['n'],
						      self.ebq_global['dx_f'],
						      self.starJacobian)
    def getConservationFluxPWC(self):
        self.e['conservationCorrectionPWC'].flat[:]=0.0
        self.getConservationResidualPWC()
        if self.updateConservationJacobian:
            self.getConservationJacobianPWC()
            self.updateConservationJacobian = False
        self.conservationSolver.solve(u=self.e['conservationCorrectionPWC'],b=self.e['conservationResidual'])
        self.e['conservationCorrectionPWC']*=-1.0
        self.updateElementBoundaryCoefficients()
        self.getConservationResidualPWC()
    def getConservationFluxPWL(self):
        self.starU.flat[:]=0.0
        self.getConservationResidualPWL()
        if self.updateConservationJacobian:
	    self.getConservationJacobianPWL()
	    self.updateConservationJacobian = False
        femIntegrals.calculateConservationFluxPWL(self.mesh.nNodes_global,
                                                  self.nNodes_internal,
                                                  self.mesh.nElements_node,
                                                  self.nodeStarOffset,
                                                  self.nodeStarJacobianOffset,
                                                  self.internalNodesArray,
                                                  self.starR,
                                                  self.starJacobian,
                                                  self.starU)
        self.getConservationResidualPWL()
    def writeBoundaryTermsEnsight(self,filename):
        caseOut=open(filename+'.case','a')
        caseOut.write('measured: '+filename+'_elementBoundaryQuadrature.geo\n')
        caseOut.write('VARIABLE\n')
        caseOut.write('vector per measured node: unitNormal '+filename+'_unitNormal.vec\n')
        caseOut.write('vector per measured node: velocityAverage '+filename+'_velocityAverage.vec\n')
        caseOut.write('vector per measured node: cnv '+filename+'_cnv.vec\n')
        caseOut.write('scalar per measured node: velocityJump '+filename+'_velocityJump.vec\n')
        caseOut.write('scalar per measured node: fluxAverage '+filename+'_fluxAverage.vec\n')
        caseOut.write('scalar per measured node: conservationFluxPWC '+filename+'_conservationFluxPWC.vec\n')
        caseOut.close()
        quadratureOut=open(filename+'_elementBoundaryQuadrature.geo','w')
        quadratureOut.write('quadrature points\n'+'particle coordinates\n')
        quadratureOut.write('%8i\n' % (self.mesh.nElementBoundaries_global*self.nElementBoundaryQuadraturePoints_elementBoundary,))
        normalOut = open(filename+'_unitNormal.vec','w')
        normalOut.write('element boundary unit normal\n')
        velocityAverageOut = open(filename+'_velocityAverage.vec','w')
        velocityAverageOut.write('element boundary velocityAverage\n')
        cnvOut = open(filename+'_cnv.vec','w')
        cnvOut.write('element boundary conservative normal velocity\n')
        velocityJumpOut = open(filename+'_velocityJump.vec','w')
        velocityJumpOut.write('element boundary velocityJump\n')
        fluxAverageOut = open(filename+'_fluxAverage.vec','w')
        fluxAverageOut.write('element boundary fluxAverage\n')
        conservationFluxPWCOut = open(filename+'_conservationFluxPWC.vec','w')
        conservationFluxPWCOut.write('element boundary conservationFluxPWC\n')
        nColumnsV = 0
        nColumnsS = 0
        for ebN in range(self.mesh.nElementBoundaries_global):
            for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    qN = 1+ ebN*self.nElementBoundaryQuadraturePoints_elementBoundary + k
                    quadratureOut.write('%8i%12.5e%12.5e%12.5e\n' % (qN,
                                                                     self.ebq_global['x'][ebN,k,0],
                                                                     self.ebq_global['x'][ebN,k,1],
                                                                     self.ebq_global['x'][ebN,k,2]))
                    if self.nSpace_global == 1:
                        normalOut.write('%12.5e%12.5e%12.5e' %(self.ebq_global['n'][ebN,k,0],
                                                               0.0,
                                                               0.0))
                        velocityAverageOut.write('%12.5e%12.5e%12.5e' %(self.ebq_global['velocityAverage'][ebN,k,0],
                                                                        0.0,
                                                                        0.0))
                        cnvOut.write('%12.5e%12.5e%12.5e' %(self.eb_global['conservationFluxPWC'][ebN]*self.ebq_global['n'][ebN,k,0],
                                                            0.0,
                                                            0.0))
                    elif self.nSpace_global == 2:
                        normalOut.write('%12.5e%12.5e%12.5e' %(self.ebq_global['n'][ebN,k,0],
                                                               self.ebq_global['n'][ebN,k,1],
                                                               0.0))
                        velocityAverageOut.write('%12.5e%12.5e%12.5e' %(self.ebq_global['velocityAverage'][ebN,k,0],
                                                                        self.ebq_global['velocityAverage'][ebN,k,1],
                                                                        0.0))
                        cnvOut.write('%12.5e%12.5e%12.5e' %(self.eb_global['conservationFluxPWC'][ebN]*self.ebq_global['n'][ebN,k,0],
                                                            self.eb_global['conservationFluxPWC'][ebN]*self.ebq_global['n'][ebN,k,1],
                                                            0.0))
                    elif self.nSpace_global == 3:
                        normalOut.write('%12.5e%12.5e%12.5e' %(self.ebq_global['n'][ebN,k,0],
                                                               self.ebq_global['n'][ebN,k,1],
                                                               self.ebq_global['n'][ebN,k,2]))
                        velocityAverageOut.write('%12.5e%12.5e%12.5e' %(self.ebq_global['velocityAverage'][ebN,k,0],
                                                                        self.ebq_global['velocityAverage'][ebN,k,1],
                                                                        self.ebq_global['velocityAverage'][ebN,k,2]))
                        cnvOut.write('%12.5e%12.5e%12.5e' %(self.eb_global['conservationFluxPWC'][ebN]*self.ebq_global['n'][ebN,k,0],
                                                            self.eb_global['conservationFluxPWC'][ebN]*self.ebq_global['n'][ebN,k,1],
                                                            self.eb_global['conservationFluxPWC'][ebN]*self.ebq_global['n'][ebN,k,2]))
                    nColumnsV+=3
                    if nColumnsV == 6:
                        nColumnsV = 0
                        normalOut.write('\n')
                        velocityAverageOut.write('\n')
                        cnvOut.write('\n')
                    velocityJumpOut.write('%12.5e' % (self.ebq_global['velocityJump'][ebN,k],))
                    #just put the average of the element boundary flux at all the quadrature points
                    fluxAverageOut.write('%12.5e' % (self.eb_global['fluxAverage'][ebN],))
                    conservationFluxPWCOut.write('%12.5e' % (self.eb_global['conservationFluxPWC'][ebN],))
                    nColumnsS+=1
                    if nColumnsS == 6:
                        nColumnsS = 0
                        velocityJumpOut.write('\n')
                        fluxAverageOut.write('\n')
                        conservationFluxPWCOut.write('\n')
        normalOut.write('\n')
        velocityAverageOut.write('\n')
        cnvOut.write('\n')
        velocityJumpOut.write('\n')
        fluxAverageOut.write('\n')
        conservationFluxPWCOut.write('\n')
        normalOut.close()
        quadratureOut.close()
        velocityAverageOut.close()
        cnvOut.close()
        velocityJumpOut.close()
        fluxAverageOut.close()
        conservationFluxPWCOut.close()
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
        if self.numericalFlux != None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                for ii in range(self.l2g['nFreeDOF'][left_eN_global]):
                    i = self.l2g['freeLocal'][left_eN_global,ii]
                    left_I = self.l2g['freeGlobal'][left_eN_global,ii]
                    if not columnIndecesDict.has_key(left_I):
                        columnIndecesDict[left_I]=set()
                    for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
                        j = self.l2g['freeLocal'][left_eN_global,jj]
                        left_J = self.l2g['freeGlobal'][left_eN_global,jj]
                        columnIndecesDict[left_I].add(left_J)
                        #jacobian[left_I,left_J]   =1.0
                    for jj in range(self.l2g['nFreeDOF'][right_eN_global]):
                        j = self.l2g['freeLocal'][right_eN_global,jj]
                        right_J = self.l2g['freeGlobal'][right_eN_global,jj]
                        columnIndecesDict[left_I].add(right_J)
                        #jacobian[left_I,right_J]  =1.0
                for ii in range(self.l2g['nFreeDOF'][right_eN_global]):
                    i = self.l2g['freeLocal'][right_eN_global,ii]
                    right_I = self.l2g['freeGlobal'][right_eN_global,ii]
                    if not columnIndecesDict.has_key(right_I):
                        columnIndecesDict[right_I]=set()
                    for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
                        j = self.l2g['freeLocal'][left_eN_global,jj]
                        left_J = self.l2g['freeGlobal'][left_eN_global,jj]
                        columnIndecesDict[right_I].add(left_J)
                        #jacobian[right_I,left_J]  =1.0
                    for  jj  in range(self.l2g['nFreeDOF'][right_eN_global]):
                        j = self.l2g['freeLocal'][right_eN_global,jj]
                        right_J = self.l2g['freeGlobal'][right_eN_global,jj]
                        columnIndecesDict[right_I].add(right_J)
                        #jacobian[right_I,right_J] =1.0
        if self.numericalFlux != None  or  (self.fluxBoundaryConditions == 'outFlow'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                for  ii  in range(self.l2g['nFreeDOF'][eN_global]):
                    i = self.l2g['freeLocal'][eN_global,ii]
                    I = self.l2g['freeGlobal'][eN_global,ii]
                    if not columnIndecesDict.has_key(I):
                        columnIndecesDict[I]=set()
                    for  jj  in range(self.l2g['nFreeDOF'][eN_global]):
                        j = self.l2g['freeLocal'][eN_global,jj]
                        J = self.l2g['freeGlobal'][eN_global,jj]
                        columnIndecesDict[I].add(J)
                        #jacobian[I,J] =1.0
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
        if self.numericalFlux != None:
            for ebNI in range(self.mesh.nInteriorElementBoundaries_global):
                ebN = self.mesh.interiorElementBoundariesArray[ebNI]
                left_eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                right_eN_global  = self.mesh.elementBoundaryElementsArray[ebN,1]
                left_ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                right_ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,1]
                for ii in range(self.l2g['nFreeDOF'][left_eN_global]):
                    i = self.l2g['freeLocal'][left_eN_global,ii]
                    left_I = self.l2g['freeGlobal'][left_eN_global,ii]
                    for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
                        j = self.l2g['freeLocal'][left_eN_global,jj]
                        left_J = self.l2g['freeGlobal'][left_eN_global,jj]
                        self.csrColumnOffsets_eb[ebN,0,0,ii,jj] = columnOffsetDict[(left_I,left_J)]
                    for jj in range(self.l2g['nFreeDOF'][right_eN_global]):
                        j = self.l2g['freeLocal'][right_eN_global,jj]
                        right_J = self.l2g['freeGlobal'][right_eN_global,jj]
                        self.csrColumnOffsets_eb[ebN,0,1,ii,jj] = columnOffsetDict[(left_I,right_J)]
                for ii in range(self.l2g['nFreeDOF'][right_eN_global]):
                    i = self.l2g['freeLocal'][right_eN_global,ii]
                    right_I = self.l2g['freeGlobal'][right_eN_global,ii]
                    for jj in range(self.l2g['nFreeDOF'][left_eN_global]):
                        j = self.l2g['freeLocal'][left_eN_global,jj]
                        left_J = self.l2g['freeGlobal'][left_eN_global,jj]
                        self.csrColumnOffsets_eb[ebN,1,0,ii,jj] = columnOffsetDict[(right_I,left_J)]
                    for  jj  in range(self.l2g['nFreeDOF'][right_eN_global]):
                        j = self.l2g['freeLocal'][right_eN_global,jj]
                        right_J = self.l2g['freeGlobal'][right_eN_global,jj]
                        self.csrColumnOffsets_eb[ebN,1,1,ii,jj] = columnOffsetDict[(right_I,right_J)]
        if self.numericalFlux != None  or (self.fluxBoundaryConditions == 'outFlow'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
                ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
                for  ii  in range(self.l2g['nFreeDOF'][eN_global]):
                    i = self.l2g['freeLocal'][eN_global,ii]
                    I = self.l2g['freeGlobal'][eN_global,ii]
                    for  jj  in range(self.l2g['nFreeDOF'][eN_global]):
                        j = self.l2g['freeLocal'][eN_global,jj]
                        J = self.l2g['freeGlobal'][eN_global,jj]
                        self.csrColumnOffsets_eb[ebN,0,0,ii,jj] = columnOffsetDict[(I,J)]
        self.nNonzerosInJacobian = len(columnOffsetDict)
        #print self.csrRowIndeces
        #print self.csrColumnOffsets
        #done building csr data structures
        return self.jacobian
    #mwf add for uniform interface with OneLevelVector
    def saveSolution(self):
        pass

class MultilevelScalarTransport:
    """Nonlinear ADR on a multilevel mesh"""
    
    def __init__(self,nd,mlMesh,
                 TrialSpaceType,
                 TestSpaceType,
                 matType,
                 dirichletConditionSetter,
                 coefficients,
                 quadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditions='noFlow',#'outFlow'
                 stabilization='2',
                 shockCapturing='1',
                 shockCapturingDiffusion=0.1,
                 conservativeFlux=None,
                 numericalFlux=None,
                 TimeIntegrationClass = BackwardEuler, #mwf add tOrder
                 tOrder = 1):

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
        print "Building ScalarTransport for each mesh"
        for mesh in mlMesh.meshList:
            print "Generating Trial Space"
            trialSpace = TrialSpaceType(mesh,nd)
            self.trialSpaceList.append(trialSpace)
            print "Generating Test Space"
            testSpace = TestSpaceType(mesh,nd)
            print "Allocating u"
            u = FiniteElementFunction(trialSpace)
            print "Allocating phi"
            phi = FiniteElementFunction(trialSpace)
            print "Setting Boundary Conditions"
            dirichletConditions=DOFBoundaryConditions(
                trialSpace,dirichletConditionSetter)
            self.bcList.append(dirichletConditions)
            print "Initializing OneLevelScalarTransport"
            scalarTransport=OneLevelScalarTransport(u,
                                                    phi,
                                                    testSpace,
                                                    matType,
                                                    dirichletConditions,
                                                    copy.deepcopy(coefficients),
                                                    quadrature,
                                                    elementBoundaryQuadrature,
                                                    fluxBoundaryConditions,
                                                    stabilization,
                                                    shockCapturing,
                                                    shockCapturingDiffusion,
                                                    conservativeFlux,
                                                    numericalFlux,
                                                    TimeIntegrationClass,
                                                    tOrder)
            self.modelList.append(scalarTransport)
            print "Allocating Vectors"
            self.u = Numeric.zeros((scalarTransport.dim,),Numeric.Float)
            self.du = Numeric.zeros((scalarTransport.dim,),Numeric.Float)
            self.r = Numeric.zeros((scalarTransport.dim,),Numeric.Float)
            self.matType = matType
            print "Allocating Jacobian"
            self.jacobian = scalarTransport.initializeJacobian()
            self.jacobianList.append(self.jacobian)
            self.uList.append(self.u)
            self.duList.append(self.du)
            self.rList.append(self.r)
        print "Building Mesh Transfers"
        self. meshTransfers = MultilevelProjectionOperators(
            mlMesh,
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
            m.inflowBoundaryBC[:]=m.inflowBoundary
    #end updateTimeHistory
    #mwf added
    def updateStage(self):
        for m,u,r in zip(self.modelList,self.uList,self.rList):
            m.timeIntegration.updateStage()

## @}

if __name__ == '__main__':
    import profile
    import pstats
    import ScalarTransportTests
    profile.run('ScalarTransportTests.runTests()','ScalarTransportProf')
    p = pstats.Stats('ScalarTransportProf')
    p.sort_stats('cumulative').print_stats(20)
    p.sort_stats('time').print_stats(20)
