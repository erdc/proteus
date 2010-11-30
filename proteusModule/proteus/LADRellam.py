from Transport import *
#import cLADR1Dellam,cLADR2Dellam,ctracking,Quadrature
import cellam,ctracking,Quadrature

"""
TODO
  high:
  med:
    need better way to set ellam specific options?
  low:
    add coupling between components 
    decide if should get rid of numerical flux use or keep?
  1D 

  2D

  3D
  
  General
     
     look at extra unknowns at outflow boundary for large time step
      and debug small oscillations
     worth trying something to get smoother inflow boundary approx. for
     step function  (e.g. add strategic time integration points)?
   
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
        #mwf debug
        #import pdb
        #pdb.set_trace()
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

        ###for tracking
        #quadrature points
        self.x_track = {}; self.t_track={}; self.t_depart={}; self.dt_track={}; self.flag_track={}; self.element_track={}
        for ci in range(self.nc):
            self.x_track[ci] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
            self.t_track[ci]   = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.t_depart[ci]   = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.dt_track[ci]   = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
            self.flag_track[ci]   = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'i')
            self.flag_track[ci].fill(-1)
            self.element_track[ci]   = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'i')
        #interpolation points for solution values
        self.x_track_ip = {}; self.t_track_ip={}; self.t_depart_ip={}; self.u_track_ip={}; self.flag_track_ip={}; self.element_track_ip={}
        self.u_dof_lim = {}; self.u_dof_lim_tmp = {} ; self.u_dof_lim_last = {}
        for ci in range(self.nc):
            self.x_track_ip[ci] = numpy.zeros((self.mesh.nElements_global,self.n_phi_ip_element[ci],3),'d')
            self.t_track_ip[ci]   = numpy.zeros((self.mesh.nElements_global,self.n_phi_ip_element[ci]),'d')
            self.t_depart_ip[ci]   = numpy.zeros((self.mesh.nElements_global,self.n_phi_ip_element[ci]),'d')
            self.u_track_ip[ci]   = numpy.zeros((self.mesh.nElements_global,self.n_phi_ip_element[ci]),'d')
            self.flag_track_ip[ci]   = numpy.zeros((self.mesh.nElements_global,self.n_phi_ip_element[ci]),'i')
            self.flag_track_ip[ci].fill(-1)
            self.element_track_ip[ci]   = numpy.zeros((self.mesh.nElements_global,self.n_phi_ip_element[ci]),'i')
            self.u_dof_lim[ci] = numpy.copy(self.u[ci].dof)
            self.u_dof_lim_tmp[ci] = numpy.copy(self.u[ci].dof)
            self.u_dof_lim_last[ci] = numpy.copy(self.u[ci].dof)
            
        self.elementBoundaryOuterNormalsArray = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,self.nSpace_global),'d')
        #todo clean this up, add set from options for zeroTol etc
        assert 'particleTracking' in dir(options), "ELLAM requires particleTracking type to be set in n file"
        self.particleTrackingType = options.particleTracking
        self.particle_tracker = options.particleTracking(self.mesh,self.nSpace_global,
                                                         activeComponentList=range(self.nc))
                       
        self.particle_tracker.setFromOptions(options)
        self.particle_tracker.updateTransportInformation(self)
        self.zeroSolutionTol_track = {}
        if 'zeroSolutionTol_track' in dir(options):
            for ci in range(self.nc):
                self.zeroSolutionTol_track[ci]=options.zeroSolutionTol_track[ci]
        else:
            for ci in range(self.nc):
                self.zeroSolutionTol_track[ci]=1.0e-8
            
        self.useBackwardTrackingForOldMass = False
        if 'useBackwardTrackingForOldMass' in dir(options):
            self.useBackwardTrackingForOldMass = options.useBackwardTrackingForOldMass

        if self.useBackwardTrackingForOldMass:
            #need this to evaluate coefficients at tracked points
            self.q_track = {}
            #deep_keys = [('u',0),('m',0),('dm',0,0),('f',0),('df',0,0),('a',0,0),('velocity',0)]
            deep_keys = set([('u',ci) for ci in range(self.nc)])
            deep_keys |= set([('m',ci) for ci in self.coefficients.mass.keys()])
            deep_keys |= set([('f',ci) for ci in self.coefficients.advection.keys()])
            deep_keys |= set([('velocity',ci) for ci in range(self.nc)])
            for ci,cjDict in self.coefficients.mass.iteritems():
                deep_keys |= set([('dm',ci,cj) for cj in cjDict.keys()])
            for ci,cjDict in self.coefficients.advection.iteritems():
                deep_keys |= set([('df',ci,cj) for cj in cjDict.keys()])
            for ci,ckDict in self.coefficients.diffusion.iteritems():
                deep_keys |= set([('a',ci,ck) for ck in ckDict.keys()])
                
            shallow_keys=set(['x'])
            for k in self.q.keys():
                if k in deep_keys:
                    self.q_track[k] = numpy.copy(self.q[k])
                elif k in shallow_keys:
                    self.q_track[k] = self.q[k]
            self.u_dof_track = {}
            for ci in range(self.nc):
                self.u_dof_track[ci] = numpy.copy(self.u[ci].dof)    
        else:
            self.q_track=self.q #could only grab shallow copies of specific keys
            self.u_dof_track = {}
            for ci in range(self.nc):
                self.u_dof_track[ci] = self.u[ci].dof
        #boundary point tracking (inflow approx)
        self.ebqe_x_track = {}; self.ebqe_t_track={}; self.ebqe_t_depart={};  self.ebqe_flag_track={}; self.ebqe_element_track={}
        for ci in range(self.nc):
            self.ebqe_x_track[ci] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
            self.ebqe_t_track[ci]   = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.ebqe_t_depart[ci]   = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.ebqe_flag_track[ci]   = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
            self.ebqe_flag_track[ci].fill(-1)
            self.ebqe_element_track[ci]   = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')

        #keep track of when to update tracking
        self.needToTrackPoints = True;
        self.tForLastTrackingStep = 0.0
        #outflow boundary approximation via trapezoidal rule
        for ci in range(self.nc):
            self.ebqe[('outflow_flux',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.ebqe[('outflow_flux_last',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        
            self.ebqe[('u',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
            self.ebqe[('grad(u)',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')


            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
            self.ebqe[('advectiveFlux_bc',ci)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        #
        self.needEBQ = options.needEBQ #could need for analytical velocity evaluation with RT0,BDM

        #data structures for slumping
        self.slumpingFlag =0
        self.rightHandSideForLimiting = {}; self.elementResidualTmp = {}; self.elementModifiedMassMatrixCorrection = {}; self.elementSlumpingParameter = {}
        for ci in range(self.nc):
            self.rightHandSideForLimiting[ci]= numpy.zeros((self.nFreeDOF_global[ci],),'d')
            self.elementResidualTmp[ci] = numpy.zeros((self.mesh.nElements_global,self.nDOF_test_element[ci]),'d')
            self.elementModifiedMassMatrixCorrection[ci] = numpy.zeros((self.mesh.nElements_global,self.nDOF_test_element[ci],self.nDOF_test_element[ci]),'d')
            self.elementSlumpingParameter[ci] = numpy.zeros((self.mesh.nElements_global,),'d')
        #
        if 'slumpingFlag' in dir(options):
            self.slumpingFlag = options.slumpingFlag
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


    def calculateElementResidualOriginalWorks(self):
        """Calculate standard element residuals needed for ellam approximation to linear ADR example"""
        import pdb
        for ci in range(self.nc):
            self.elementResidual[ci].fill(0.0)
        #calculate normal Eulerian weak terms
        # (m^{n+1,w^{n+1}) + (\Delta t(x) a\grad u, grad w^{n+1}) + (\Delta t(x) r,w^{n+1}) 
        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck in ckDict.keys():
                #weight by time step size
                self.q[('dt*grad(w)*dV_a',ck,ci)][:] = self.q[('grad(w)*dV_a',ck,ci)]
                #todo need faster loop
                for j in range(self.nDOF_trial_element[0]):
                    for I in range(self.nSpace_global):
                        self.q[('dt*grad(w)*dV_a',ck,ci)][:,:,j,I] *= self.dt_track[ci]
                if self.sd:
                    cfemIntegrals.updateDiffusion_weak_sd(self.coefficients.sdInfo[(ci,ck)][0],self.coefficients.sdInfo[(ci,ck)][1],
                                                          self.q[('a',ci,ck)],
                                                          self.q[('grad(phi)',ck)],
                                                          self.q[('dt*grad(w)*dV_a',ck,ci)],
                                                          self.elementResidual[ci])
                else:
                    cfemIntegrals.updateDiffusion_weak_lowmem(self.q[('a',ci,ck)],
                                                              self.q[('grad(phi)',ck)],
                                                              self.q[('dt*grad(w)*dV_a',ck,ci)],
                                                              self.elementResidual[ci])
                        

        for ci in self.coefficients.reaction.keys():
            #weight by time step size
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #print "LADRellam evalres rxn= %s " % (self.q[('r',0)])
            self.q[('dt*w*dV_r',ci)][:] = self.q[('w*dV_r',ci)]
            #todo need faster loop
            for j in range(self.nDOF_trial_element[0]):
                self.q[('dt*w*dV_r',ci)][:,:,j]   *= self.dt_track[ci]
            cfemIntegrals.updateReaction_weak(self.q[('r',ci)],
                                              self.q[('dt*w*dV_r',ci)],
                                              self.elementResidual[ci])
        for ci in self.coefficients.mass.keys():
            #note not dm/dt but just m
            cfemIntegrals.updateMass_weak(self.q[('m',ci)],
                                          self.q[('w*dV_m',ci)],
                                          self.elementResidual[ci])

        #
        # (m^{n},w^{n+1})
        self.approximateOldMassIntegral(self.elementResidual)
        #inflow
        self.approximateInflowBoundaryIntegral(self.elementResidual)
        #outflow
        self.approximateOutflowBoundaryIntegral()
    def calculateElementResidual(self):
        """
        Calculate standard element residuals needed for ellam approximation to linear ADR example
        Switch order around to facilitate slumping/limiting, compute explicit/rhs parts first
        """
        import pdb
       
        for ci in range(self.nc):
            self.elementResidual[ci].fill(0.0)
        if self.slumpingFlag == 1:
            for ci in range(self.nc):
                self.elementResidualTmp[ci].fill(0.0)
        for ci in self.coefficients.reaction.keys():
            #weight by time step size
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #print "LADRellam evalres rxn= %s " % (self.q[('r',0)])
            self.q[('dt*w*dV_r',ci)][:] = self.q[('w*dV_r',ci)]
            #todo need faster loop
            for j in range(self.nDOF_trial_element[0]):
                self.q[('dt*w*dV_r',ci)][:,:,j]   *= self.dt_track[ci]
            cfemIntegrals.updateReaction_weak(self.q[('r',ci)],
                                              self.q[('dt*w*dV_r',ci)],
                                              self.elementResidual[ci])
        
        #
        # (m^{n},w^{n+1})
        self.approximateOldMassIntegral(self.elementResidual)
        #inflow
        self.approximateInflowBoundaryIntegral(self.elementResidual)
        #outflow
        self.approximateOutflowBoundaryIntegral()
        if self.slumpingFlag == 1:
            for ci in range(self.nc):
                self.elementResidualTmp[ci] -= self.elementResidual[ci]
        
        # (m^{n+1,w^{n+1}) + (\Delta t(x) a\grad u, grad w^{n+1}) + (\Delta t(x) r,w^{n+1}) 
        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck in ckDict.keys():
                #weight by time step size
                self.q[('dt*grad(w)*dV_a',ck,ci)][:] = self.q[('grad(w)*dV_a',ck,ci)]
                #todo need faster loop
                for j in range(self.nDOF_trial_element[0]):
                    for I in range(self.nSpace_global):
                        self.q[('dt*grad(w)*dV_a',ck,ci)][:,:,j,I] *= self.dt_track[ci]
                if self.sd:
                    cfemIntegrals.updateDiffusion_weak_sd(self.coefficients.sdInfo[(ci,ck)][0],self.coefficients.sdInfo[(ci,ck)][1],
                                                          self.q[('a',ci,ck)],
                                                          self.q[('grad(phi)',ck)],
                                                          self.q[('dt*grad(w)*dV_a',ck,ci)],
                                                          self.elementResidual[ci])
                else:
                    cfemIntegrals.updateDiffusion_weak_lowmem(self.q[('a',ci,ck)],
                                                              self.q[('grad(phi)',ck)],
                                                              self.q[('dt*grad(w)*dV_a',ck,ci)],
                                                              self.elementResidual[ci])
                        

        for ci in self.coefficients.mass.keys():
            #note not dm/dt but just m
            cfemIntegrals.updateMass_weak(self.q[('m',ci)],
                                          self.q[('w*dV_m',ci)],
                                          self.elementResidual[ci])

        #mwf debug
        #pdb.set_trace()
        if self.slumpingFlag == 1:
            for ci in range(self.nc):
                #assemble right hand side vector
                self.rightHandSideForLimiting[ci].fill(0.)
                cfemIntegrals.updateGlobalResidualFromElementResidual(self.offset[ci],
                                                                      self.stride[ci],
                                                                      self.l2g[ci]['nFreeDOF'],
                                                                      self.l2g[ci]['freeLocal'],
                                                                      self.l2g[ci]['freeGlobal'],
                                                                      self.elementResidualTmp[ci],
                                                                      self.rightHandSideForLimiting[ci]);
            #calculate element level lumping parameters and
            #subtract off element level mass correction from residual
            if self.nSpace_global == 1:
                cellam.calculateSlumpedMassApproximation1d(self.u[ci].femSpace.dofMap.l2g,
                                                           self.mesh.elementNeighborsArray,
                                                           self.u[ci].dof,self.u[ci].dof,
                                                           self.q[('dm',ci,ci)],
                                                           self.q[('w',ci)],
                                                           self.q[('v',ci)],
                                                           self.q[('dV_u',ci)],
                                                           self.rightHandSideForLimiting[ci],
                                                           self.elementResidual[ci],
                                                           self.elementSlumpingParameter[ci],
                                                           self.elementModifiedMassMatrixCorrection[ci])
                
            elif self.nSpace_global == 2:
                cellam.calculateSlumpedMassApproximation2d(self.u[ci].femSpace.dofMap.l2g,
                                                           self.mesh.elementNeighborsArray,
                                                           self.u[ci].dof,self.u[ci].dof,
                                                           self.q[('dm',ci,ci)],
                                                           self.q[('w',ci)],
                                                           self.q[('v',ci)],
                                                           self.q[('dV_u',ci)],
                                                           self.rightHandSideForLimiting[ci],
                                                           self.elementResidual[ci],
                                                           self.elementSlumpingParameter[ci],
                                                           self.elementModifiedMassMatrixCorrection[ci])
                
        elif self.slumpingFlag == 2:
            #start by using current solution to do limiting, then try back tracking
            if self.nSpace_global == 1:
                #mwf debug
                #import pdb
                #pdb.set_trace()
                cellam.calculateBerzinsSlumpedMassApproximation1d(self.u[ci].femSpace.dofMap.l2g,
                                                                  self.mesh.elementNeighborsArray,
                                                                  self.u[ci].dof,self.u_dof_lim[ci],
                                                                  self.q[('dm',ci,ci)],
                                                                  self.q[('w',ci)],
                                                                  self.q[('v',ci)],
                                                                  self.q[('dV_u',ci)],
                                                                  self.rightHandSideForLimiting[ci],
                                                                  self.elementResidual[ci],
                                                                  self.elementModifiedMassMatrixCorrection[ci])
                
    def calculateElementJacobian(self):
        for ci in range(self.nc):
            for cj in self.coefficients.stencil[ci]:
                self.elementJacobian[ci][cj].fill(0.0)
        ##\todo optimize nonlinear diffusion Jacobian calculation for the  different combinations of nonlinear a and phi
        
        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                for cj in set(cjDict.keys()+self.coefficients.potential[ck].keys()):
                    #assume dt weighting has been set already
                    if self.sd:
                        cfemIntegrals.updateDiffusionJacobian_weak_sd(self.coefficients.sdInfo[(ci,ck)][0],self.coefficients.sdInfo[(ci,ck)][1],
                                                                      self.phi[ck].femSpace.dofMap.l2g,
                                                                      self.q[('a',ci,ck)],
                                                                      self.q[('da',ci,ck,cj)],
                                                                      self.q[('grad(phi)',ck)],
                                                                      self.q[('dt*grad(w)*dV_a',ck,ci)],
                                                                      self.dphi[(ck,cj)].dof,
                                                                      self.q[('v',cj)],
                                                                      self.q[('grad(v)',cj)],
                                                                      self.elementJacobian[ci][cj])
                    else:
                        cfemIntegrals.updateDiffusionJacobian_weak_lowmem(self.phi[ck].femSpace.dofMap.l2g,
                                                                          self.q[('a',ci,ck)],
                                                                          self.q[('da',ci,ck,cj)],
                                                                          self.q[('grad(phi)',ck)],
                                                                          self.q[('dt*grad(w)*dV_a',ck,ci)],
                                                                          self.dphi[(ck,cj)].dof,
                                                                          self.q[('v',cj)],
                                                                          self.q[('grad(v)',cj)],
                                                                          self.elementJacobian[ci][cj])
        for ci,cjDict in self.coefficients.reaction.iteritems():
            for cj in cjDict:
                #assume dt weighting has been set already
                cfemIntegrals.updateReactionJacobian_weak_lowmem(self.q[('dr',ci,cj)],
                                                                 self.q[('v',cj)],
                                                                 self.q[('dt*w*dV_r',ci)],
                                                                 self.elementJacobian[ci][cj])
                    
        for ci,cjDict in self.coefficients.mass.iteritems():
            for cj in cjDict:
                cfemIntegrals.updateMassJacobian_weak_lowmem(self.q[('dm',ci,cj)],
                                                             self.q[('v',cj)],
                                                             self.q[('w*dV_m',ci)],
                                                             self.elementJacobian[ci][cj])
                
        if self.slumpingFlag == 1:
            #mwf debug
            #import pdb
            #pdb.set_trace()
            useC = True
            for ci,cjDict in self.coefficients.mass.iteritems():
                for cj in cjDict:
                    if useC:
                        cellam.updateElementJacobianWithSlumpedMassApproximation(self.elementSlumpingParameter[ci],
                                                                                 self.elementJacobian[ci][cj])
                    else:
                        for eN in range(self.mesh.nElements_global):
                            for i in range(self.nDOF_test_element[ci]):
                                self.elementJacobian[ci][cj][eN,i,i] += (self.nDOF_trial_element[cj]-1)*self.elementSlumpingParameter[ci][eN]
                                for j in range(i):
                                    self.elementJacobian[ci][cj][eN,i,j] -= self.elementSlumpingParameter[ci][eN]
                                for j in range(i+1,self.nDOF_trial_element[cj]):
                                    self.elementJacobian[ci][cj][eN,i,j] -= self.elementSlumpingParameter[ci][eN]

        elif self.slumpingFlag == 2:
            #mwf debug
            #import pdb
            #pdb.set_trace()
            useC = False
            for ci,cjDict in self.coefficients.mass.iteritems():
                for cj in cjDict:
                    if useC:
                        cellam.updateElementJacobianWithSlumpedMassCorrection(self.elementModifiedMassMatrixCorrection[ci],
                                                                              self.elementJacobian[ci][cj])
                    else:
                        for eN in range(self.mesh.nElements_global):
                            for i in range(self.nDOF_test_element[ci]):
                                for j in range(self.nDOF_trial_element[cj]):
                                    self.elementJacobian[ci][cj][eN,i,j] += self.elementModifiedMassMatrixCorrection[ci][eN,i,j]

    def calculateExteriorElementBoundaryJacobian(self):
        for jDict in self.fluxJacobian_exterior.values():
            for j in jDict.values():
                j.fill(0.0)
        self.approximateOutflowBoundaryIntegralJacobian()
    def updateTimeHistory(self,T,resetFromDOF=False):
        """
        todo find a better place to make sure know when a step is done
        because if step failes need to retrack
        """
        OneLevelTransport.updateTimeHistory(self,T,resetFromDOF)
        self.needToTrackPoints = True
        for ci in range(self.nc):
            self.ebqe[('outflow_flux_last',ci)].flat[:] = self.ebqe[('outflow_flux',ci)].flat
        #todo put this in time integration
        #mwf debug
        #import pdb
        #pdb.set_trace()
        if self.useBackwardTrackingForOldMass:
            for ci in range(self.nc):
                self.u_dof_track[ci].flat[:] = self.u[ci].dof.flat
        if self.slumpingFlag == 2:
            for ci in range(self.nc):
                self.u_dof_lim_last[ci].flat[:] = self.u[ci].dof.flat
            
        log("ELLAM t= %s Global conservation= %s " % (self.timeIntegration.t,numpy.sum(self.elementResidual[0].flat)),level=0)
    def setInitialConditions(self,getInitialConditionsDict,T=0.0):
        OneLevelTransport.setInitialConditions(self,getInitialConditionsDict,T=T)
        if self.useBackwardTrackingForOldMass:
            for ci in range(self.nc):
                self.u_dof_track[0].flat[:] = self.u[0].dof.flat

    def trackQuadraturePoints(self):
        """
        track quadrature points in q['x'] backward from t^{n+1} --> t^{n}, 
          loads
             x_track[0]      : location of point at end of tracking
             t_track[0]      : time tracking ended
             flag_track[0]   : -1  -- point in interior at tOut
			       -2  -- point exited domain somewhere in (tIn,tOut)
                               -3  -- did not track (e.g., v = 0 or u = 0)
             element_track[0]     : element containing point at end of tracking
        save time steps for domain in 
             dt_track[0] = t^{n+1} - t_track[0]
        Then
        track quadrature points in q['x'] forward from t^n --> t^{n+1}, 
          save 
             q[('x_track',0)]     : location of point at end of tracking
             t_track[0]      : time tracking ended
             flag_track[0]   : -1  -- point in interior at tOut
			       -2  -- point exited domain somewhere in (tIn,tOut)
                               -3  -- did not track (e.g., v = 0 or u = 0)
             element_track[0]     : element containing point at end of tracking
             
        """
        import pdb
        #backward
        if (self.needToTrackPoints and self.timeIntegration.t > self.timeIntegration.tLast + 1.0e-8 or
            abs(self.tForLastTrackingStep-self.timeIntegration.t) > 1.0e-8):
            #mwf debug
            #pdb.set_trace()
            #update velocity fields for particle tracking
            x_depart = {}
            nPoints_track  = {}
            for ci in range(self.nc):
                self.particle_tracker.setTrackingVelocity(self.coefficients.adjoint_velocity_dofs_last[ci],ci,
                                                          self.coefficients.adjoint_velocity_times_last[ci],
                                                          timeLevel=0,
                                                          trackingVelocity_l2g=self.coefficients.adjoint_velocity_l2g[ci])
                self.particle_tracker.setTrackingVelocity(self.coefficients.adjoint_velocity_dofs[ci],ci,
                                                          self.coefficients.adjoint_velocity_times[ci],
                                                          timeLevel=1)


                log(" LADRellam tracking integration points backward ci=%s" % ci,level=2) 
                self.t_depart[ci].fill(self.timeIntegration.t)
                #in desired output time, out actual time
                self.t_track[ci].fill(self.timeIntegration.tLast)
                #try all points, now set to -1 to try, -3 to skip, 0 or greater if a node of the mesh
                self.flag_track[ci].fill(-1)
                #don't skip points with small concentration when going backwards in case near a point with nonzero grad?
                skipPointsWithZeroSolution = 0
                
                for k in range(self.element_track[ci].shape[1]):
                    self.element_track[ci][:,k] = numpy.arange(self.mesh.nElements_global,dtype='i')

                x_depart[ci] = self.q['x']
                nPoints_track[ci] = self.mesh.nElements_global*self.nQuadraturePoints_element
                if skipPointsWithZeroSolution:
                    cellam.tagNegligibleIntegrationPoints(nPoints_track[ci],
                                                          self.zeroSolutionTol_track[ci],
                                                          x_depart[ci],
                                                          self.q[('u',ci)],
                                                          self.flag_track[ci])
                    
            #todo make sure activeComponents set explicitly?
            
            self.particle_tracker.backwardTrack(self.t_depart,
                                                self.t_track,
                                                nPoints_track,
                                                x_depart,
                                                self.element_track,
                                                self.x_track,
                                                self.flag_track)


            #mwf debug
            #pdb.set_trace()
            for ci in range(self.nc):
                self.dt_track[ci]  = numpy.copy(self.t_depart[ci])
                self.dt_track[ci] -= self.t_track[ci]
                #mwf debug tracked locations
                #for eN in range(self.mesh.nElements_global):
                #    for k in range(self.nQuadraturePoints_element):
                #        eN_track = self.element_track[0][eN,k]
                #        x0=self.mesh.nodeArray[self.mesh.elementNodesArray[eN_track,0],0]
                #        x1=self.mesh.nodeArray[self.mesh.elementNodesArray[eN_track,1],0]
                #        if (self.x_track[0][eN,k][0] < min(x0,x1)-1.0e-8 or 
                #            self.x_track[0][eN,k][0] > max(x0,x1)+1.0e-8):
                #            pdb.set_trace()
                        
            if not self.useBackwardTrackingForOldMass:
                for ci in range(self.nc):
                    log(" LADRellam tracking integration points forward ci=%s " % ci,level=2) 
                    #forward
                    self.t_depart[ci].fill(self.timeIntegration.tLast)
                    self.t_track[ci].fill(self.timeIntegration.t)
                    #todo this doesn't have any effect now
                    #do skip points with small concentration when going forwards
                    skipPointsWithZeroSolution = 1
                    #try all points, now set to -1 to try, -3 to skip, 0 or greater if a node of the mesh
                    self.flag_track[ci].fill(-1)
                    #todo get rid of this loop
                    #for eN in range(self.mesh.nElements_global):
                    #    self.element_track[0][eN,:] = eN
                    #self.element_track[0].flat[:] = numpy.arange(self.mesh.nElements_global,dtype='i').repeat(self.element_track[0].shape[1],axis=None)
                    for k in range(self.element_track[ci].shape[1]):
                        self.element_track[ci][:,k] = numpy.arange(self.mesh.nElements_global,dtype='i')
                    #todo set zero tol in initialization
                    if skipPointsWithZeroSolution:
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()
                        cellam.tagNegligibleIntegrationPoints(nPoints_track[ci],
                                                              self.zeroSolutionTol_track[ci],
                                                              x_depart[ci],
                                                              self.q[('u',ci)],
                                                              self.flag_track[ci])
                #todo make sure activeComponents set explicitly?
                self.particle_tracker.forwardTrack(self.t_depart,
                                                   self.t_track,
                                                   nPoints_track,
                                                   x_depart,
                                                   self.element_track,
                                                   self.x_track,
                                                   self.flag_track)


                    #mwf debug tracked locations
                    #for eN in range(self.mesh.nElements_global):
                    #    for k in range(self.nQuadraturePoints_element):
                    #        eN_track = self.element_track[0][eN,k]
                    #        x0=self.mesh.nodeArray[self.mesh.elementNodesArray[eN_track,0],0]
                    #        x1=self.mesh.nodeArray[self.mesh.elementNodesArray[eN_track,1],0]
                    #        if (self.x_track[0][eN,k][0] < min(x0,x1)-1.0e-8 or 
                    #            self.x_track[0][eN,k][0] > max(x0,x1)+1.0e-8):
                    #            pdb.set_trace()

            needToTrackInterpolationPoints = (self.slumpingFlag == 2)
            if needToTrackInterpolationPoints:
                x_depart_ip = {}
                nPoints_track_ip  = {}
                for ci in range(self.nc):
                    #todo switch this to characteristic velocity, need _last values!
                    self.particle_tracker.setTrackingVelocity(self.coefficients.adjoint_velocity_dofs_last[ci],ci,
                                                              self.coefficients.adjoint_velocity_times_last[ci],
                                                              timeLevel=0,
                                                              trackingVelocity_l2g=self.coefficients.adjoint_velocity_l2g[ci])
                    self.particle_tracker.setTrackingVelocity(self.coefficients.adjoint_velocity_dofs[ci],ci,
                                                              self.coefficients.adjoint_velocity_times[ci],
                                                              timeLevel=1)


                    log(" LADRellam tracking integration points backward ci=%s" % ci,level=2) 
                    self.t_depart_ip[ci].fill(self.timeIntegration.t)
                    #in desired output time, out actual time
                    self.t_track_ip[ci].fill(self.timeIntegration.tLast)
                    #try all points, now set to -1 to try, -3 to skip, 0 or greater if a node of the mesh
                    self.flag_track_ip[ci].fill(-1)
                    #don't skip points with small concentration when going backwards in case near a point with nonzero grad?
                    skipPointsWithZeroSolution = 0

                    for k in range(self.element_track_ip[ci].shape[1]):
                        self.element_track_ip[ci][:,k] = numpy.arange(self.mesh.nElements_global,dtype='i')

                    x_depart_ip[ci] = self.u[ci].femSpace.interpolationPoints
                    nPoints_track_ip[ci] = self.mesh.nElements_global*x_depart_ip[ci].shape[1]
                    if skipPointsWithZeroSolution:
                        assert False, "todo need u array that's the same shape as interpolation points array"
                        #cellam.tagNegligibleIntegrationPoints(nPoints_track_ip[ci],
                        #                                      self.zeroSolutionTol_track[ci],
                        #                                      x_depart_ip[ci],
                        #                                      self.q[('u',ci)],
                        #                                      self.flag_track_ip[ci])

                #todo make sure activeComponents set explicitly?

                self.particle_tracker.backwardTrack(self.t_depart_ip,
                                                    self.t_track_ip,
                                                    nPoints_track_ip,
                                                    x_depart_ip,
                                                    self.element_track_ip,
                                                    self.x_track_ip,
                                                    self.flag_track_ip)

                for ci in range(self.nc):
                    #evaluate solution at tracked interpolation points using old time level dofs
                    cellam.evaluateSolutionAtTrackedPoints(self.nSpace_global,
                                                           self.nDOF_trial_element[ci],
                                                           nPoints_track_ip[ci],
                                                           self.mesh.nElements_global,
                                                           self.mesh.nNodes_global,
                                                           self.mesh.nNodes_element,
                                                           self.mesh.nElementBoundaries_element,
                                                           self.mesh.nodeArray,
                                                           self.mesh.elementNodesArray,
                                                           self.mesh.elementNeighborsArray,
                                                           self.elementBoundaryOuterNormalsArray,
                                                           self.x_track_ip[ci],
                                                           self.t_track_ip[ci],
                                                           self.element_track_ip[ci],
                                                           self.flag_track_ip[ci],
                                                           self.u[ci].femSpace.dofMap.l2g,
                                                           self.u_dof_lim_last[ci],#todo put this in time integration?
                                                           self.u_track_ip[ci])
                    #use finite element machinery to project solution at new time level
                    self.u_dof_lim_tmp[ci][:] = self.u[ci].dof
                    self.u[ci].projectFromInterpolationConditions(self.u_track_ip[ci])
                    self.u_dof_lim[ci][:] = self.u[ci].dof
                    self.u[ci].dof[:] = self.u_dof_lim_tmp[ci]
                
            #end tracking interpolation points
            self.needToTrackPoints = False
            self.tForLastTrackingStep=self.timeIntegration.t
            #mwf debug
            #pdb.set_trace()
        #

    def approximateOldMassIntegral(self,elementRes):
        """
        approximate weak integral
        \int_{\Omega} m^{n} w^{n+1} \dV using forward tracking
        """
        if self.useBackwardTrackingForOldMass:
            log("LADRellam evaluating old mass integral with backtracking",level=2)
            #assumes that x_track, t_track etc correctly set from backtracking step in trackQuadraturePoints
            for ci in range(self.nc):
                cellam.evaluateSolutionAtTrackedPoints(self.nSpace_global,
                                                      self.nDOF_trial_element[ci],
                                                      self.mesh.nElements_global*self.nQuadraturePoints_element,
                                                      self.mesh.nElements_global,
                                                      self.mesh.nNodes_global,
                                                      self.mesh.nNodes_element,
                                                      self.mesh.nElementBoundaries_element,
                                                      self.mesh.nodeArray,
                                                      self.mesh.elementNodesArray,
                                                      self.mesh.elementNeighborsArray,
                                                      self.elementBoundaryOuterNormalsArray,
                                                      self.x_track[ci],
                                                      self.t_track[ci],
                                                      self.element_track[ci],
                                                      self.flag_track[ci],
                                                      self.u[ci].femSpace.dofMap.l2g,
                                                      self.u_dof_track[ci],#todo put this in time integration?
                                                      self.q_track[('u',ci)])

            #now evaluate as a standard mass integral
            #todo get rid of all of this, just want mass 
            self.q_track['dV']= self.q['dV']
            #mwf debug
            #import pdb
            #pdb.set_trace()
            self.coefficients.evaluate(self.timeIntegration.tLast,self.q_track)
            for ci in range(self.nc):
                #have to scale by -1
                self.q_track[('m',ci)] *= -1.
                cfemIntegrals.updateMass_weak(self.q_track[('m',ci)],
                                              self.q[('w*dV_m',ci)],
                                              elementRes[ci])
                self.q_track[('m',ci)] *= -1.
            #mwf debug
            #import pdb
            #pdb.set_trace()
        else:
            log("LADRellam evaluating old mass integral with forwardtracking",level=2)
            #mwf debug
            #import pdb
            #pdb.set_trace()
            for ci in range(self.nc):
                cellam.updateOldMass_weak(self.nSpace_global,
                                          self.nDOF_test_element[ci],
                                          self.mesh.nElements_global,
                                          self.mesh.nNodes_global,
                                          self.mesh.nNodes_element,
                                          self.mesh.nElementBoundaries_element,
                                          self.nQuadraturePoints_element,
                                          self.mesh.nodeArray,
                                          self.mesh.elementNodesArray,
                                          self.mesh.elementNeighborsArray,
                                          self.elementBoundaryOuterNormalsArray,
                                          self.q['dV'],
                                          self.x_track[ci],
                                          self.t_track[ci],
                                          self.element_track[ci],
                                          self.flag_track[ci],
                                          self.u[ci].femSpace.dofMap.l2g,
                                          self.timeIntegration.m_last[ci],
                                          elementRes[ci])

    def approximateInflowBoundaryIntegral(self,elementRes):
        """
        approximate term

         \int_{t^n}^{t^{n+1}}  \int_{\Gamma_{I}\sigma^b w \dS \dt
        
        numerically using composite trapezoidal rule in time (and space too)

          \sum_{p=1}^{NT}\sum_{q=1}^{N_{q,b}}\Delta t^{p}\sigma^b(x_{q},t^p)w^{n+1}_{i}(\tilde{x}_q,t^{n+1})} W_q

        Here (x_q,t^p) tracks forward to  (\tilde{x}_q,t^{n+1}) and w^{n+1}_{i} is any test function with support
          covering (\tilde{x}_q,t^{n+1})

        only points on inflow boundary are tracked  
        """
        if self.timeIntegration.t > self.timeIntegration.tLast + 1.0e-8:
            #mwf debug
            #import pdb
            #pdb.set_trace()
            #update velocity fields for particle tracking
            ebqe_x_depart = {}
            ebqe_nPoints_track  = {}
            for ci in range(self.nc):
                self.particle_tracker.setTrackingVelocity(self.coefficients.adjoint_velocity_dofs_last[ci],ci,
                                                          self.coefficients.adjoint_velocity_times_last[ci],
                                                          timeLevel=0,
                                                          trackingVelocity_l2g=self.coefficients.adjoint_velocity_l2g[ci])
                self.particle_tracker.setTrackingVelocity(self.coefficients.adjoint_velocity_dofs[ci],ci,
                                                          self.coefficients.adjoint_velocity_times[ci],
                                                          timeLevel=1)
                ebqe_nPoints_track[ci]=self.mesh.nExteriorElementBoundaries_global*self.nElementBoundaryQuadraturePoints_elementBoundary
                ebqe_x_depart[ci] = self.ebqe['x']
            self.NT = max(2,4*int(ceil(self.timeIntegration.runCFL)))
            dtp = (self.timeIntegration.t-self.timeIntegration.tLast)/float(self.NT)
            integrationTimes = numpy.arange(self.NT+1,dtype='d')*dtp + self.timeIntegration.tLast
            integrationTimeWeights=numpy.zeros(self.NT+1,'d'); integrationTimeWeights.fill(dtp)
            integrationTimeWeights[0] *= 0.5; integrationTimeWeights[-1] *= 0.5

            for tpi,dtpi in zip(integrationTimes,integrationTimeWeights):
                for ci in range(self.nc):
                    #figure out which points on inflow need to be tracked
                    cellam.markInflowBoundaryPoints(self.nSpace_global,
                                                    self.timeIntegration.tLast,
                                                    self.timeIntegration.t,
                                                    tpi,
                                                    self.mesh.nExteriorElementBoundaries_global,
                                                    self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                    self.mesh.exteriorElementBoundariesArray,
                                                    self.mesh.elementBoundaryElementsArray,
                                                    self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                    self.ebqe['x'],
                                                    self.ebqe['n'],
                                                    self.coefficients.ebqe[('velocity',ci)],#need to have time varying v
                                                    self.coefficients.ebqe[('velocity',ci)],
                                                    self.numericalFlux.isDOFBoundary[ci],
                                                    self.ebqe[('advectiveFlux_bc_flag',ci)],
                                                    self.ebqe_element_track[ci],
                                                    self.ebqe_flag_track[ci])
                    
                    #track forward
                    self.ebqe_t_depart[ci].fill(tpi)
                    self.ebqe_t_track[ci].fill(self.timeIntegration.t)
                    #need to skip points with small boundary flux when tracking inflow boundary
                    skipPointsWithZeroSolution = 1

                    if skipPointsWithZeroSolution:
                        #mwf debug
                        #import pdb
                        #pdb.set_trace()
                        cellam.tagNegligibleIntegrationPoints(ebqe_nPoints_track[ci],
                                                              self.zeroSolutionTol_track[ci],
                                                              ebqe_x_depart[ci],
                                                              self.ebqe[('advectiveFlux_bc',ci)],
                                                              self.ebqe_flag_track[ci])
                    
                direction = 1.0 #forward tracking
                if self.timeIntegration.t > tpi + 1.0e-8: 
                    self.particle_tracker.forwardTrack(self.ebqe_t_depart,
                                                           self.ebqe_t_track,
                                                           ebqe_nPoints_track,
                                                           ebqe_x_depart,
                                                           self.ebqe_element_track,
                                                           self.ebqe_x_track,
                                                           self.ebqe_flag_track)

                    for ci in range(self.nc):
                        #accumulate into correct locations in residual
                        cellam.accumulateInflowFlux(self.nSpace_global,
                                                    self.nDOF_test_element[ci],
                                                    self.mesh.nElements_global,
                                                    self.mesh.nNodes_global,
                                                    self.mesh.nNodes_element,
                                                    self.mesh.nElementBoundaries_element,
                                                    self.mesh.nExteriorElementBoundaries_global,
                                                    self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                    self.mesh.nodeArray,
                                                    self.mesh.elementNodesArray,
                                                    self.mesh.elementNeighborsArray,
                                                    self.mesh.exteriorElementBoundariesArray,
                                                    self.mesh.elementBoundaryElementsArray,
                                                    self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                    self.elementBoundaryOuterNormalsArray,
                                                    tpi,
                                                    dtpi,
                                                    self.ebqe['dS'],
                                                    self.ebqe_x_track[ci],
                                                    self.ebqe_t_track[ci],
                                                    self.ebqe_element_track[ci],
                                                    self.ebqe_flag_track[ci],
                                                    self.u[ci].femSpace.dofMap.l2g,
                                                    self.u[ci].dof,
                                                    elementRes[ci], 
                                                    self.coefficients.sdInfo[(ci,ci)][0], #todo fix
                                                    self.coefficients.sdInfo[(ci,ci)][1],
                                                    self.ebqe[('advectiveFlux_bc_flag',ci)],
                                                    self.ebqe[('advectiveFlux_bc',ci)])


    def approximateOutflowBoundaryIntegral(self):
        """
        approximate term

         \int_{t^n}^{t^{n+1}}  \int_{\Gamma_{O}} f w \dS \dt
        
        numerically using trapezoidal rule in time 


        """
        for ci in range(self.nc):
            #accumulate into correct locations in residual
            cellam.updateExteriorOutflowBoundaryFlux(self.timeIntegration.t-self.timeIntegration.tLast,
                                                     self.nSpace_global,
                                                     self.nDOF_test_element[ci],
                                                     self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                     self.mesh.nExteriorElementBoundaries_global,
                                                     self.mesh.exteriorElementBoundariesArray,
                                                     self.mesh.elementBoundaryElementsArray,
                                                     self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                     self.coefficients.ebqe[('velocity',ci)],
                                                     self.ebqe[('n')],
                                                     self.ebqe[('outflow_flux_last',ci)],
                                                     self.ebqe[('w*dS_f',ci)],
                                                     self.ebqe[('u',ci)],
                                                     self.u[ci].femSpace.dofMap.l2g,
                                                     self.ebqe[('outflow_flux',ci)],
                                                     self.elementResidual[ci])

                        
                    
    def approximateOutflowBoundaryIntegralJacobian(self):
        """
        approximate jacobian for term

         \int_{t^n}^{t^{n+1}}  \int_{\Gamma_{O}} f w \dS \dt
        
        numerically using trapezoidal rule in time 


        """
        for ci in range(self.nc):
            cellam.updateExteriorOutflowBoundaryFluxJacobian(self.timeIntegration.t-self.timeIntegration.tLast,
                                                             self.nSpace_global,
                                                             self.nDOF_test_element[ci],
                                                             self.nDOF_trial_element[ci],
                                                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                             self.mesh.nExteriorElementBoundaries_global,
                                                             self.mesh.exteriorElementBoundariesArray,
                                                             self.mesh.elementBoundaryElementsArray,
                                                             self.mesh.elementBoundaryLocalElementBoundariesArray,
                                                             self.coefficients.ebqe[('velocity',ci)],
                                                             self.ebqe[('n')],
                                                             self.ebqe[('outflow_flux_last',ci)],
                                                             self.ebqe[('w*dS_f',ci)],
                                                             self.ebqe[('u',ci)],
                                                             self.ebqe[('v',ci)],
                                                             self.fluxJacobian_exterior[ci][ci])


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
        #this could just be a call to evaluate velocity for now
        self.calculateCoefficients()
        
        #setup points and time steps for tracking
        self.trackQuadraturePoints()

        self.calculateElementResidual()
        for ci in range(self.nc):
            cfemIntegrals.updateGlobalResidualFromElementResidual(self.offset[ci],
                                                                  self.stride[ci],
                                                                  self.l2g[ci]['nFreeDOF'],
                                                                  self.l2g[ci]['freeLocal'],
                                                                  self.l2g[ci]['freeGlobal'],
                                                                  self.elementResidual[ci],
                                                                  r);
        #handle inflow separately for now
        #self.approximateInflowBoundaryIntegral(u,r)
        #and outflow
        #self.approximateOutflowBoundaryIntegral(u,r)
        #mwf debug
        #pdb.set_trace()
        #mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1
    def getJacobian(self,jacobian):
        import superluWrappers
        import numpy
        import pdb
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
				       jacobian)

        
        #mwf debug
        #pdb.set_trace()
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
        import pdb
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
        #extra boundary normal information for 2d, 3d to save need for ebq array
        boundaryNormals = numpy.array(self.testSpace[0].elementMaps.referenceElement.boundaryUnitNormalList,dtype='d')
        ctracking.getOuterNormals_affineSimplex(boundaryNormals,
                                                self.q['inverse(J)'],
                                                self.elementBoundaryOuterNormalsArray)
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
        #mwf hack
        #TODO make sure coefficients has access to quadrature points for velocity evaluation??
        self.coefficients.elementQuadraturePoints = self.elementQuadraturePoints
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
        #mwf TODO cleanup make sure coefficients has access to quadrature points for velocity evaluation??
        self.coefficients.elementBoundaryQuadraturePoints = self.elementBoundaryQuadraturePoints
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
       
    def estimate_mt(self):
        pass

    def calculateElementBoundaryCoefficients(self):
        """
        Calculate the nonlinear coefficients at the element boundary quadrature points
        """
        pass
