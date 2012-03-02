from Transport import *
#import cLADR1Dellam,cLADR2Dellam,ctracking,Quadrature
import cellam,ctracking,Quadrature
import ELLAMtools
"""
TODO
  high:
    double check Kuzmin FCT why not exact mass cons. if aij and mij are symmetric?
      --> looks like consistent mass matrix quadrature values are varying in the 6th digit too?
    double check sign conventions in 1d, looks like Kuzmin and Moeller still assume Fij on rhs in limiting?
    double check max min and see if cutting them off at 0,1 helps over-under shoot
      --> didn't see much differcence.
    Probably really need to use low order solution as mass lumping in limiting if want to recover
     exactly the monotinicity constraints. Would Kuzmin 06 limiting be different?

    fix Kuzmin FCT implementation so that works for systems (that is the jacobian)
    check why sometimes seg faults on 1d slug test on first rerun
  med:
    need better way to set ellam specific options?
    move setupInitialElementLocations to c, allow for different types of tracking points
  low:
    add coupling between components
    decide if should get rid of numerical flux use or keep?
  1D

  2D

  3D

  General
     try some strategic integration point approximations in 1d, 2d, 3d, --> now testing
     slumping in 2d,3d --> now testing 2d ideas,
     Kuzmin Turek approach 1d,2d,3d --> now testing Kuzmin_Moeller_etal_10 approach

     look at extra unknowns at outflow boundary for large time step
      and debug small oscillations
     worth trying something to get smoother inflow boundary approx. for
     step function  (e.g. add strategic time integration points)?


For SSIPs type approach ...

   see if picking SSIPs at new time level as a refined quadrature approximation for elements with
     stuff going on but still defining quadrature at old time level helps
   Do we need to go to a full refinement scheme?

"""
class OneLevelLADR(OneLevelTransport):
    nCalls=0
    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='defaultName',
                 reuse_trial_and_test_quadrature=True,
                 sd = True,
                 movingDomain=False):
        """
        Allocate storage and initialize some variables.

        uDict   -- a dictionary of FiniteElementFunction objects

        phiDict -- a dictionary of FiniteElementFunction objects

        testSpaceDict -- a dictionary of FiniteElementSpace objects

        dofBoundaryConditionsDict -- a dictionary of DOFBoundaryConditions objects for
        the Dirichlet conditions

        coefficients -- a TransportCoefficients object

        elementQuadratureDict -- a dictionary of dictionaries of quadrature rules for each
        element integral in each component equation

        elementBoundaryQuadratureDict -- a dictionary of dictionaries of quadrature rules
        for each element boundary integral in each component equation

        stabilization

        shockCapturing

        numericalFlux

        The constructor sets the input arguments, calculates
        dimensions, and allocates storage. The meanings of variable
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
        q          -- at element quadrature, unique to elements
        ebq        -- at element boundary quadrature, unique to elements
        ebq_global -- at element boundary quadrature, unique to element boundary
        ebqe       -- at element boundary quadrature, unique to global, exterior element boundary
        phi_ip     -- at the generalized interpolation points required to build a nonlinear  phi
        """
        #
        #set the objects describing the method and boundary conditions
        #
        self.movingDomain=movingDomain
        self.tLast_mesh=None
        #
        self.name=name
        self.sd=sd
        self.Hess=False
        self.lowmem=True
        self.timeTerm=True#allow turning off  the  time derivative
        #self.lowmem=False
        self.testIsTrial=True
        self.phiTrialIsTrial=True
        self.u = uDict
        self.ua = {}#analytical solutions
        self.phi  = phiDict
        self.dphi={}
        for ck,phi in phiDict.iteritems():
            if coefficients.potential.has_key(ck):
                for cj in coefficients.potential[ck].keys():
                    self.dphi[(ck,cj)] = FiniteElementFunction(phi.femSpace)
            else:
                self.dphi[(ck,ck)] = FiniteElementFunction(phi.femSpace)
        #check for nonlinearities in the diffusion coefficient that don't match the potential
        for ci,ckDict in coefficients.diffusion.iteritems():
            #for ck,cjDict in coefficients.diffusion.iteritems(): #cek: bug?
            for ck,cjDict in ckDict.iteritems():
                for cj in cjDict.keys():
                    if not self.dphi.has_key((ck,cj)):
                        self.dphi[(ck,cj)] = FiniteElementFunction(phi.femSpace)
        self.matType = matType
        #try to reuse test and trial information across components if spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature#True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1,coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        ## Simplicial Mesh
        self.mesh = self.u[0].femSpace.mesh #assume the same mesh for  all components for now
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        self.dirichletNodeSetList=None #explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.conservativeFlux = conservativeFluxDict #no velocity post-processing for now
        self.fluxBoundaryConditions=fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict=advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        #determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        #cek come back
        if self.stabilization != None:
            for ci in range(self.nc):
                if coefficients.mass.has_key(ci):
                    for flag in coefficients.mass[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if  coefficients.advection.has_key(ci):
                    for  flag  in coefficients.advection[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if  coefficients.diffusion.has_key(ci):
                    for diffusionDict in coefficients.diffusion[ci].values():
                        for  flag  in diffusionDict.values():
                            if flag != 'constant':
                                self.stabilizationIsNonlinear=True
                if  coefficients.potential.has_key(ci):
                    for flag in coefficients.potential[ci].values():
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if coefficients.reaction.has_key(ci):
                    for flag in coefficients.reaction[ci].values():
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if coefficients.hamiltonian.has_key(ci):
                    for flag in coefficients.hamiltonian[ci].values():
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
        #determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci  in range(self.nc):
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux != None) or
                                                 (numericalFluxType != None) or
                                                 (self.fluxBoundaryConditions[ci] == 'outFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'mixedFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        #calculate some dimensions
        #
        self.nSpace_global    = self.u[0].femSpace.nSpace_global #assume same space dim for all variables
        self.nDOF_trial_element     = [u_j.femSpace.max_nDOF_element for  u_j in self.u.values()]
        self.nDOF_phi_trial_element     = [phi_k.femSpace.max_nDOF_element for  phi_k in self.phi.values()]
        self.n_phi_ip_element = [phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for  phi_k in self.phi.values()]
        self.nDOF_test_element     = [femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global  = [dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
        self.nVDOF_element    = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global)
        #
        NonlinearEquation.__init__(self,self.nFreeVDOF_global)
        #
        #build the quadrature point dictionaries from the input (this
        #is just for convenience so that the input doesn't have to be
        #complete)
        #
        elementQuadratureDict={}
        elemQuadIsDict = isinstance(elementQuadrature,dict)
        if elemQuadIsDict: #set terms manually
            for I in self.coefficients.elementIntegralKeys:
                if elementQuadrature.has_key(I):
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization != None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature
        if self.shockCapturing != None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff',ci,ci)):
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature[('numDiff',ci,ci)]
                    else:
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        elementBoundaryQuadratureDict={}
        if isinstance(elementBoundaryQuadrature,dict): #set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        #mwf include tag telling me which indices are which quadrature rule?
        (self.elementQuadraturePoints,self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element*self.mesh.nElements_global
        #
        #Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints,
         self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[0]
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global*
                                                        self.mesh.nElementBoundaries_element*
                                                        self.nElementBoundaryQuadraturePoints_elementBoundary)

        #
        #storage dictionaries
        self.scalars_element = set()
        #
        #simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q={}
        self.ebq={}
        self.ebq_global={}
        self.ebqe={}
        self.phi_ip={}
        #mesh
        self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.q['det(J)'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['abs(det(J))'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['J'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global,self.nSpace_global),'d')
        self.q['inverse(J)'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global,self.nSpace_global),'d')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.ebqe['g'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                       self.nElementBoundaryQuadraturePoints_elementBoundary,
                                       max(1,self.nSpace_global-1),
                                       max(1,self.nSpace_global-1)),
                                      'd')
        self.ebqe['inverse(J)'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global,self.nSpace_global),'d')
        self.ebqe['hat(x)'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.ebqe['bar(x)'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.ebqe['sqrt(det(g))'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('n')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #shape
        self.q[('v',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0]),'d')
        self.q[('w',0)] = self.q[('v',0)]
        self.q[('grad(v)',0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0],self.nSpace_global),'d')
        self.q[('grad(w)',0)] =  self.q[('grad(v)',0)]
        self.q[('grad(w)*dV',0)]   =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0],self.nSpace_global),'d')
        self.q[('grad(w)*dV_f',0)] = self.q[('grad(w)*dV',0)]
        #todo get rid of dV_{f,a}, etc
        self.q[('w*dV',0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0]),'d')
        self.q[('w*dV_m',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0]),'d')
        #assume all components are the same space for now
        shapeKeysForAlias = ['v','w','grad(v)','grad(w)*dV','grad(w)*dV_f','w*dV','w*dV_m']
        for ci in range(1,self.nc):
            for key in shapeKeysForAlias:
                key_ci = (key,ci)
                key_0  = (key,0)
                self.q[key_ci] = self.q[key_0]
        #ELLAM weights stiffness, body integrals by dt
        for ci in range(self.nc):
            self.q[('dt*grad(w)*dV',ci)]= numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[ci],self.nSpace_global),'d')
        #
        self.ebqe[('v',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[0]),'d')
        self.ebqe[('w',0)] = self.ebqe[('v',0)]
        self.ebqe[('grad(v)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[0],self.nSpace_global),'d')
        self.ebqe[('w*dS_f',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[0]),'d')
        #assume all components are the same space for now
        shapeKeysForAlias = ['v','w','grad(v)','w*dS_f']
        for ci in range(1,self.nc):
            for key in shapeKeysForAlias:
                key_ci = (key,ci)
                key_0  = (key,0)
                self.ebqe[key_ci] = self.ebqe[key_0]

        for ci in range(self.nc):
            self.q[('u',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.q[('grad(u)',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        #f
        for ci in self.coefficients.advection.keys():
            self.q[('f',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
            for cj in self.coefficients.advection[ci].keys():
                self.q[('df',ci,cj)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
            self.ebqe[('f',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
            for cj in self.coefficients.advection[ci].keys():
                self.ebqe[('df',ci,cj)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')

        #a, linear dispersion single component

        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                for flag in cjDict.values():
                    assert flag == 'constant', "Error potential %s LADRellam does not handle diffusion = %s yet" % (ck,flag)

                if self.coefficients.sdInfo != None and (ci,ck) in self.coefficients.sdInfo.keys():
                    self.q[('a',ci,ck)] = numpy.zeros(
                        (self.mesh.nElements_global,
                         self.nQuadraturePoints_element,
                         self.coefficients.sdInfo[(ci,ck)][0][self.nSpace_global]),
                        'd')
                    for cj in cjDict.keys():
                        self.q[('da',ci,ck,cj)] = numpy.zeros(
                            (self.mesh.nElements_global,
                             self.nQuadraturePoints_element,
                             self.coefficients.sdInfo[(ci,ck)][0][self.nSpace_global]),
                            'd')
                    self.ebqe[('a',ci,ck)]=numpy.zeros(
                        (self.mesh.nExteriorElementBoundaries_global,
                         self.nElementBoundaryQuadraturePoints_elementBoundary,
                         self.coefficients.sdInfo[(ci,ck)][0][self.nSpace_global]),
                        'd')
                    for cj in cjDict.keys():
                        self.ebqe[('da',ci,ck,cj)]=numpy.zeros(
                            (self.mesh.nExteriorElementBoundaries_global,
                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                             self.coefficients.sdInfo[(ci,ck)][0][self.nSpace_global]),
                            'd')

                else:
                    self.q[('a',ci,ck)]=numpy.zeros(
                        (self.mesh.nElements_global,
                         self.nQuadraturePoints_element,
                         self.nSpace_global,
                         self.nSpace_global),
                        'd')
                    for cj in cjDict.keys():
                        self.q[('da',ci,ck,cj)]=numpy.zeros(
                            (self.mesh.nElements_global,
                             self.nQuadraturePoints_element,
                             self.nSpace_global,
                             self.nSpace_global),
                            'd')
                    self.ebqe[('a',ci,ck)]=numpy.zeros(
                        (self.mesh.nExteriorElementBoundaries_global,
                         self.nElementBoundaryQuadraturePoints_elementBoundary,
                         self.nSpace_global,
                         self.nSpace_global),
                        'd')
                    for cj in cjDict.keys():
                        self.ebqe[('da',ci,ck,cj)]=numpy.zeros(
                            (self.mesh.nExteriorElementBoundaries_global,
                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                             self.nSpace_global,
                             self.nSpace_global),
                            'd')
                #dense storage
            self.q[('grad(w)*dV_a',ci,ck)]   = self.q[('grad(w)*dV_f',ci)]
            self.q[('dt*grad(w)*dV_a',ci,ck)]= self.q[('dt*grad(w)*dV',ci)]
        #ci,ckDict
        #linear potential only for now, need to change for e.g., Buckley Leverett
        for ck in self.phi.keys():
            self.phi[ck].dof[:]=self.u[ck].dof
            self.q[('grad(phi)',ck)] = self.q[('grad(u)',ck)]
        for key in self.dphi.keys():
            self.dphi[key].dof.fill(1.0)
            self.q[('dphi',key[0],key[1])] = numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')



#         if self.coefficients.diffusion.has_key(0):
#             for ck,flag in self.coefficients.diffusion[0][0].iteritems():
#                 assert  self.coefficients.diffusion[0][0][ck] == 'constant', "Error potential %s LADRellam does not handle diffusion = %s yet" % (ck,flag)
#             if self.coefficients.sdInfo != None and (0,0) in self.coefficients.sdInfo.keys():
#                 self.q[('a',0,0)] = numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.nQuadraturePoints_element,
#                      self.coefficients.sdInfo[(0,0)][0][self.nSpace_global]),
#                     'd')
#                 self.q[('da',0,0,0)] = numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.nQuadraturePoints_element,
#                      self.coefficients.sdInfo[(0,0)][0][self.nSpace_global]),
#                     'd')
#                 self.ebqe[('a',0,0)]=numpy.zeros(
#                     (self.mesh.nExteriorElementBoundaries_global,
#                      self.nElementBoundaryQuadraturePoints_elementBoundary,
#                      self.coefficients.sdInfo[(0,0)][0][self.nSpace_global]),
#                     'd')
#                 self.ebqe[('da',0,0,0)]=numpy.zeros(
#                     (self.mesh.nExteriorElementBoundaries_global,
#                      self.nElementBoundaryQuadraturePoints_elementBoundary,
#                      self.coefficients.sdInfo[(0,0)][0][self.nSpace_global]),
#                     'd')

#             else:
#                 self.q[('a',0,0)]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.nQuadraturePoints_element,
#                      self.nSpace_global,
#                      self.nSpace_global),
#                     'd')
#                 self.q[('da',0,0,0)]=numpy.zeros(
#                     (self.mesh.nElements_global,
#                      self.nQuadraturePoints_element,
#                      self.nSpace_global,
#                      self.nSpace_global),
#                     'd')
#                 self.ebqe[('a',0,0)]=numpy.zeros(
#                     (self.mesh.nExteriorElementBoundaries_global,
#                      self.nElementBoundaryQuadraturePoints_elementBoundary,
#                      self.nSpace_global,
#                      self.nSpace_global),
#                     'd')
#                 self.ebqe[('da',0,0,0)]=numpy.zeros(
#                     (self.mesh.nExteriorElementBoundaries_global,
#                      self.nElementBoundaryQuadraturePoints_elementBoundary,
#                      self.nSpace_global,
#                      self.nSpace_global),
#                     'd')
#             #
#             self.phi[0].dof[:]=self.u[0].dof
#             self.dphi[(0,0)].dof.fill(1.0)
#             self.q[('grad(phi)',0)] = self.q[('grad(u)',0)]
#             self.q[('dphi',0,0)] = numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')

#             self.q[('grad(w)*dV_a',0,0)]   = self.q[('grad(w)*dV_f',0)]
#             self.q[('dt*grad(w)*dV_a',0,0)]= self.q[('dt*grad(w)*dV',0)]

        #r 'constant' ie not a function of solution but go ahead and include dr for now
        for ci,cjDict in self.coefficients.reaction.iteritems():
            self.q[('r',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            for cj in cjDict.keys():
                self.q[('dr',ci,cj)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.q[('w*dV_r',ci)] = self.q[('w*dV',ci)]
            self.q[('dt*w*dV_r',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0]),'d')
            self.ebqe[('r',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')

        #m
        for ci,cjDict in self.coefficients.mass.iteritems():
            self.q[('m',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            for cj in cjDict.keys():
                self.q[('dm',ci,cj)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.q[('mt',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.q[('m_last',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.q[('m_tmp',ci)] = self.q[('m',ci)]
            self.q[('cfl',ci)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.q[('numDiff',ci,ci)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.ebqe[('m',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            for cj in cjDict.keys():
                self.ebqe[('dm',ci,cj)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')


        ###ellam specific options with defauls here
        self.ellamDiscretization = ELLAMtools.ELLAMdiscretization(self,options)

        #
        self.needEBQ = options.needEBQ #could need for analytical velocity evaluation with RT0,BDM

        #beg normal stuff allocating things
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()

        if self.needEBQ:
            for k in ['x','hat(x)']:
                self.ebq[k] = numpy.zeros((self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_element,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary,
                                           3),'d')
            self.ebq['n'] = numpy.zeros((self.mesh.nElements_global,
                                         self.mesh.nElementBoundaries_element,
                                         self.nElementBoundaryQuadraturePoints_elementBoundary,
                                         self.nSpace_global),'d')
            self.ebq['inverse(J)'] = numpy.zeros((self.mesh.nElements_global,
                                                  self.mesh.nElementBoundaries_element,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                  self.nSpace_global,
                                                  self.nSpace_global),'d')
            #allocate the metric tensor
            self.ebq['g'] = numpy.zeros((self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_element,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary,
                                           max(1,self.nSpace_global-1),
                                           max(1,self.nSpace_global-1)),
                                          'd')
            log(memory("element boundary quadrature","LADRellam"),level=4)
            ebq_keys = ['sqrt(det(g))']
            ebq_keys.extend([('u',ci) for ci in range(self.nc)])
            for k in ebq_keys:
                self.ebq[k] = numpy.zeros((self.mesh.nElements_global,
                                           self.mesh.nElementBoundaries_element,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary),'d')

            #test and trial info
            self.ebq[('w',0)] = numpy.zeros((self.mesh.nElements_global,
                                             self.mesh.nElementBoundaries_element,
                                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                                             self.nDOF_trial_element[0]),'d')
            for ci in range(1,self.nc):
                self.ebq[('w',ci)] = self.ebq[('w',0)]
            for ci in range(self.nc):
                self.ebq[('v',ci)] = self.ebq[('w',0)]

            #ebq_global info
            self.ebq_global['x'] = numpy.zeros((self.mesh.nElementBoundaries_global,
                                                self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                3),'d')
            self.ebq_global['n'] = numpy.zeros((self.mesh.nElementBoundaries_global,
                                                self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                self.nSpace_global),'d')
        #
        # allocate residual and Jacobian storage
        #
        self.elementResidual = [numpy.zeros(
            (self.mesh.nElements_global,
             self.nDOF_test_element[ci]),
            'd') for ci in range(self.nc)]
        self.elementSpatialResidual = [numpy.zeros(
            (self.mesh.nElements_global,
             self.nDOF_test_element[ci]),
            'd') for ci in range(self.nc)]
        self.elementJacobian = {}
        for ci in range(self.nc):
            self.elementJacobian[ci]={}
            for cj in range(self.nc):
                if cj in self.coefficients.stencil[ci]:
                    self.elementJacobian[ci][cj] = numpy.zeros(
                        (self.mesh.nElements_global,
                         self.nDOF_test_element[ci],
                         self.nDOF_trial_element[cj]),
                        'd')
        #
        self.fluxJacobian_exterior = {}
        for ci in range(self.nc):
            self.fluxJacobian_exterior[ci]={}
            for cj in self.coefficients.stencil[ci]:
                self.fluxJacobian_exterior[ci][cj] = numpy.zeros(
                    (self.mesh.nExteriorElementBoundaries_global,
                     self.nElementBoundaryQuadraturePoints_elementBoundary,
                     self.nDOF_trial_element[cj]),
                    'd')

        #
        #
        #
        #
        log(memory("element and element boundary Jacobians","OneLevelTransport"),level=4)
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
            self.inflowBoundaryBC[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,),'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nDOF_trial_element[cj]),'d')
            self.inflowFlux[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        #identify the internal nodes this is ought to be in mesh
        ##\todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global   = self.mesh.elementBoundaryElementsArray[ebN,0]
            ebN_element  = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN,0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global,i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray=numpy.zeros((self.nNodes_internal,),'i')
        for nI,n in enumerate(self.internalNodes):
            self.internalNodesArray[nI]=n
        #
        del self.internalNodes
        self.internalNodes = None
        log("Updating local to global mappings",2)
        self.updateLocal2Global()
        log("Building time integration object",2)
        log(memory("inflowBC, internalNodes,updateLocal2Global","OneLevelTransport"),level=4)
        #mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(self,integrateInterpolationPoints=True)
        else:
            self.timeIntegration = TimeIntegrationClass(self)

        if options != None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration","OneLevelTransport"),level=4)
        log("Calculating numerical quadrature formulas",2)
        self.calculateQuadrature()
        #lay out components/equations contiguously for now
        self.offset = [0]
        for ci in range(1,self.nc):
            self.offset += [self.offset[ci-1]+self.nFreeDOF_global[ci-1]]
        self.stride = [1 for ci in range(self.nc)]
        #use contiguous layout of components for parallel, requires weak DBC's
        comm = Comm.get()
        self.comm=comm
        if comm.size() > 1:
            assert numericalFluxType != None and numericalFluxType.useWeakDirichletConditions,"You must use a numerical flux to apply weak boundary conditions for parallel runs"
            self.offset = [0]
            for ci in range(1,self.nc):
                self.offset += [ci]
            self.stride = [self.nc for ci in range(self.nc)]
        #
        log(memory("stride+offset","OneLevelTransport"),level=4)
        if numericalFluxType != None:
            if options == None or options.periodicDirichletConditions == None:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict,
                                                       options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        #set penalty terms
        #cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN,k] = self.numericalFlux.penalty_constant/(self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        #penalty term
        #cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE,k] = self.numericalFlux.penalty_constant/self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux","OneLevelTransport"),level=4)
        self.elementEffectiveDiametersArray  = self.mesh.elementInnerDiametersArray
        #use post processing tools to get conservative fluxes, None by default
        import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)
        log(memory("velocity postprocessor","OneLevelTransport"),level=4)
        #helper for writing out data storage
        import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        #TODO get rid of this
        for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc',ci)].shape,'i')
            for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1

        if hasattr(self.numericalFlux,'setDirichletValues'):
            self.numericalFlux.setDirichletValues(self.ebqe)
        if not hasattr(self.numericalFlux,'isDOFBoundary'):
            self.numericalFlux.isDOFBoundary = {}
            for ci in range(self.nc):
                self.numericalFlux.isDOFBoundary[ci]= numpy.zeros(self.ebqe[('u',ci)].shape,'i')
        if not hasattr(self.numericalFlux,'ebqe'):
            self.numericalFlux.ebqe = {}
            for ci in range(self.nc):
                self.numericalFlux.ebqe[('u',ci)]= numpy.zeros(self.ebqe[('u',ci)].shape,'d')

    def calculateElementCoefficients(self):
        """
        calculate the nonlinear coefficients at the quadrature points and nodes
        this version is simplified to eliminate unnecessary logic, can always recover
        the full version by using base class
        """
        #
        #get u,grad(u), and grad(u)Xgrad(w) at the quadrature points
        #
        for cj in range(self.nc):
            self.u[cj].getValues(self.q[('v',cj)],
                                 self.q[('u',cj)])
            if self.q.has_key(('grad(u)',cj)):
                self.u[cj].getGradientValues(self.q[('grad(v)',cj)],
                                             self.q[('grad(u)',cj)])
        #
        #get functions of (t,x,u) at the quadrature points
        #
        self.coefficients.evaluate(self.timeIntegration.t,self.q)
        log("Coefficients on element",level=10,data=self.q)
        #
        # time integration is handled directly in ELLAM weak approximation, don't have a hook for
        # doing that via a time integration object (could if it were a direct Lagrange Galerkin formulation I believe)
        # however, need to set time integration's m_tmp if use that anywhere
        #if self.timeTerm:
        #    self.timeIntegration.calculateElementCoefficients(self.q)

        #todo eventually can add nonlinear potential here

        #cek and mwf need to go through this section to clean up, some of next two blocks could go to calcQuad
        #
        #todo need non-diagonal dependence?
        for ci in range(self.nc):
            cfemIntegrals.calculateCFLADR(self.elementEffectiveDiametersArray,
                                          self.q[('dm',ci,ci)],
                                          self.q[('df',ci,ci)],#could just be velocity
                                          self.q[('cfl',ci)])


    def calculateExteriorElementBoundaryCoefficients(self):
        """
        Calculate the nonlinear coefficients at global exterior element boundary quadrature points
        this version has simplified logic to reflect linear advection dispersion reaction for ellam
        """
        #
        #get u and grad(u) at the quadrature points
        #
        for ci in range(self.nc):
            self.u[ci].getValuesGlobalExteriorTrace(self.ebqe[('v',ci)],self.ebqe[('u',ci)])
            if self.ebqe.has_key(('grad(u)',ci)):
                self.u[ci].getGradientValuesGlobalExteriorTrace(self.ebqe[('grad(v)',ci)],self.ebqe[('grad(u)',ci)])
        #
        #get coefficients at the element boundary quadrature points
        #
        self.coefficients.evaluate(t = self.timeIntegration.t, c = self.ebqe)
        #
        #time integration, handled directly in ELLAM formulation
        #
        #ignore numerical flux for now
        #if self.numericalFlux != None:
        #    self.numericalFlux.calculateExteriorNumericalFlux(self.inflowFlag,self.q,self.ebqe)
        #flux boundary conditions specified through advective flux
        for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
            for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1



    def calculateElementResidual(self):
        """
        Calculate standard element residuals needed for ellam approximation to linear ADR example
        Switch order around to facilitate slumping/limiting, compute explicit/rhs parts first
        """
        import pdb

        for ci in range(self.nc):
            self.elementResidual[ci].fill(0.0)
        #
        self.ellamDiscretization.updateElementResidual(self.elementResidual)



    def calculateElementJacobian(self):
        for ci in range(self.nc):
            for cj in self.coefficients.stencil[ci]:
                self.elementJacobian[ci][cj].fill(0.0)

        self.ellamDiscretization.updateElementJacobian(self.elementJacobian)
    def calculateExteriorElementBoundaryJacobian(self):
        for jDict in self.fluxJacobian_exterior.values():
            for j in jDict.values():
                j.fill(0.0)
        self.ellamDiscretization.updateExteriorElementBoundaryJacobian(self.fluxJacobian_exterior)

    def updateTimeHistory(self,T,resetFromDOF=False):
        """
        todo find a better place to make sure know when a step is done
        because if step failes need to retrack
        """
        OneLevelTransport.updateTimeHistory(self,T,resetFromDOF)
        self.nonlinear_function_evaluations = 0
        self.ellamDiscretization.updateTimeHistory(T,resetFromDOF)
    def setInitialConditions(self,getInitialConditionsDict,T=0.0):
        OneLevelTransport.setInitialConditions(self,getInitialConditionsDict,T=T)
        self.ellamDiscretization.setInitialConditions(getInitialConditionsDict,T=T)

    def getResidual(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        #mwf debug
        #pdb.set_trace()
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        #cek can put in logic to skip of BC's don't depend on t or u
        #Dirichlet boundary conditions
        #if hasattr(self.numericalFlux,'setDirichletValues'):
        self.numericalFlux.setDirichletValues(self.ebqe)
        self.calculateCoefficients()


        self.calculateElementResidual()
        for ci in range(self.nc):
            cfemIntegrals.updateGlobalResidualFromElementResidual(self.offset[ci],
                                                                  self.stride[ci],
                                                                  self.l2g[ci]['nFreeDOF'],
                                                                  self.l2g[ci]['freeLocal'],
                                                                  self.l2g[ci]['freeGlobal'],
                                                                  self.elementResidual[ci],
                                                                  r);
        #
        self.nonlinear_function_evaluations += 1
    def getJacobian(self,jacobian):
        import superluWrappers
        import numpy
        import pdb
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)


        self.calculateElementJacobian()
        self.calculateExteriorElementBoundaryJacobian()
        for ci in range(self.nc):
            for cj in self.coefficients.stencil[ci]:
                #
                #element contributions from standard Eulerian integrals
                #
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[ci]['nFreeDOF'],
                                                                          self.l2g[ci]['freeLocal'],
                                                                          self.l2g[cj]['nFreeDOF'],
                                                                          self.l2g[cj]['freeLocal'],
                                                                          self.csrRowIndeces[(ci,cj)],
                                                                          self.csrColumnOffsets[(ci,cj)],
                                                                          self.elementJacobian[ci][cj],
                                                                          jacobian)

                cfemIntegrals.updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(self.mesh.exteriorElementBoundariesArray,
                                                                                              self.mesh.elementBoundaryElementsArray,
                                                                                              self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                                              self.l2g[ci]['nFreeDOF'],
                                                                                              self.l2g[ci]['freeLocal'],
                                                                                              self.l2g[cj]['nFreeDOF'],
                                                                                              self.l2g[cj]['freeLocal'],
                                                                                              self.csrRowIndeces[(ci,cj)],
                                                                                              self.csrColumnOffsets_eb[(ci,cj)],
                                                                                              self.fluxJacobian_exterior[ci][cj],
                                                                                              self.ebqe[('w*dS_f',ci)],
                                                                                              jacobian)
        #outflow boundary flux terms contribute too
        #self.approximateOutflowBoundaryIntegralGlobalJacobian(jacobian)


        log("Jacobian ",level=10,data=jacobian)
        #mwf debug
        #jacobian.fwrite("matdebug.txt")
        #pdb.set_trace()
        #mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian
    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        #
        #get physical locations of quadrature points and jacobian information there
        #assume all components live on the same mesh
        #
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
                                                  self.q['x'])
        if self.movingDomain:
            if self.tLast_mesh != None:
                self.q['xt'][:]=self.q['x']
                self.q['xt']-=self.q['x_last']
                alpha = 1.0/(self.t_mesh - self.tLast_mesh)
                self.q['xt']*=alpha
            else:
                self.q['xt'][:]=0.0
            self.q['x_last'][:]=self.q['x']
        self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                         self.q['J'],
                                                         self.q['inverse(J)'],
                                                         self.q['det(J)'])
        self.q['abs(det(J))']=numpy.absolute(self.q['det(J)'])
        #
        # get physical space integration weights
        #
        self.q['dV'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        cfemIntegrals.calculateIntegrationWeights(self.q['abs(det(J))'],
                                                  self.elementQuadratureWeights[('u',0)],
                                                  self.q['dV'])
        for ci in range(self.nc): self.q[('dV_u',ci)] = self.q['dV']
        #
        #get shape information at the quadrature points
        #
        self.testSpace[0].getBasisValues(self.elementQuadraturePoints,
                                         self.q[('w',0)])
        cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u',0)],
                                             self.q['abs(det(J))'],
                                             self.q[('w',0)],
                                             self.q[('w*dV',0)])
        cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('m',0)],
                                             self.q['abs(det(J))'],
                                             self.q[('w',0)],
                                             self.q[('w*dV_m',0)])
        self.testSpace[0].getBasisGradientValues(self.elementQuadraturePoints,
                                                  self.q['inverse(J)'],
                                                  self.q[('grad(w)',0)])
        cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('u',0)],
                                                      self.q['abs(det(J))'],
                                                      self.q[('grad(w)',0)],
                                                      self.q[('grad(w)*dV',0)])

        #
        self.ellamDiscretization.updateElementQuadrature(self.q)
        #
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)

    def calculateElementBoundaryQuadrature(self):
        if self.needEBQ:
            #
            #get physical locations of element boundary quadrature points
            #
            #assume all components live on the same mesh
            self.u[0].femSpace.elementMaps.getValuesTrace(self.elementBoundaryQuadraturePoints,
                                                          self.ebq['x'])

            self.u[0].femSpace.elementMaps.getJacobianValuesTrace(self.elementBoundaryQuadraturePoints,
                                                                  self.ebq['inverse(J)'],
                                                                  self.ebq['g'],
                                                                  self.ebq['sqrt(det(g))'],
                                                                  self.ebq['n'])
            useC=True
            cfemIntegrals.copyLeftElementBoundaryInfo(self.mesh.elementBoundaryElementsArray,
                                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                      self.mesh.exteriorElementBoundariesArray,
                                                      self.mesh.interiorElementBoundariesArray,
                                                      self.ebq['x'],
                                                      self.ebq['n'],
                                                      self.ebq_global['x'],
                                                      self.ebq_global['n'])
            self.u[0].femSpace.elementMaps.getInverseValuesTrace(self.ebq['inverse(J)'],self.ebq['x'],self.ebq['hat(x)'])
            self.u[0].femSpace.elementMaps.getPermutations(self.ebq['hat(x)'])

            for ci in range(self.nc):
                if self.ebq.has_key(('dS_u',ci)):
                    cfemIntegrals.calculateElementBoundaryIntegrationWeights(self.ebq['sqrt(det(g))'],
                                                                             self.elementBoundaryQuadratureWeights[('u',ci)],
                                                                             self.ebq[('dS_u',ci)])
            self.coefficients.initializeElementBoundaryQuadrature(self.timeIntegration.t,self.ebq,self.ebq_global)


    def calculateExteriorElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        #
        #get physical locations of element boundary quadrature points
        #
        #assume all components live on the same mesh
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        #
        #get metric tensor and unit normals
        #
        if self.movingDomain:
            if self.tLast_mesh != None:
                self.ebqe['xt'][:]=self.ebqe['x']
                self.ebqe['xt']-=self.ebqe['x_last']
                alpha = 1.0/(self.t_mesh - self.tLast_mesh)
                self.ebqe['xt']*=alpha
            else:
                self.ebqe['xt'][:]=0.0
            self.ebqe['x_last'][:]=self.ebqe['x']
            self.u[0].femSpace.elementMaps.getJacobianValuesGlobalExteriorTrace_movingDomain(self.elementBoundaryQuadraturePoints,
                                                                                             self.ebqe['xt'],
                                                                                             self.ebqe['inverse(J)'],
                                                                                             self.ebqe['g'],
                                                                                             self.ebqe['sqrt(det(g))'],
                                                                                             self.ebqe['n'])
        else:
            self.u[0].femSpace.elementMaps.getJacobianValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                                self.ebqe['inverse(J)'],
                                                                                self.ebqe['g'],
                                                                                self.ebqe['sqrt(det(g))'],
                                                                                self.ebqe['n'])
        #now map the physical points back to the reference element
        #assume all components live  on same mesh
        self.u[0].femSpace.elementMaps.getInverseValuesGlobalExteriorTrace(self.ebqe['inverse(J)'],self.ebqe['x'],self.ebqe['hat(x)'])
        #
        #since the points on the reference boundary may be reordered on many right element boundaries, we
        #have to use an array of reference boundary points on all element boundaries
        #first copy the left reference element boundary quadrature points from the reference element boundary
        self.testSpace[0].getBasisValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                             self.ebqe[('w',0)])
        cfemIntegrals.calculateWeightedShapeGlobalExteriorTrace(self.mesh.exteriorElementBoundariesArray,
                                                                self.mesh.elementBoundaryElementsArray,
                                                                self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                                self.elementBoundaryQuadratureWeights[('f',0)],
                                                                self.ebqe['sqrt(det(g))'],
                                                                self.ebqe[('w',0)],
                                                                self.ebqe[('w*dS_f',0)])
        self.u[0].femSpace.getBasisGradientValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                      self.ebqe['inverse(J)'],
                                                                      self.ebqe[('grad(v)',0)])
        #setup flux boundary conditions
        self.fluxBoundaryConditionsObjectsDict = dict([(cj,FluxBoundaryConditions(self.mesh,
                                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                  self.ebqe[('x')],
                                                                                  self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.ebqe['dS'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        cfemIntegrals.calculateIntegrationWeights(self.ebqe['sqrt(det(g))'],
                                                  self.elementBoundaryQuadratureWeights[('u',0)],
                                                  self.ebqe['dS'])
        for ci in range(self.nc): self.ebqe[('dS',ci)] = self.ebqe['dS']
        #
        self.ellamDiscretization.calculateExteriorElementBoundaryQuadrature(self.ebqe)
        #
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)


    def estimate_mt(self):
        pass

    def calculateElementBoundaryCoefficients(self):
        """
        Calculate the nonlinear coefficients at the element boundary quadrature points
        """
        pass
