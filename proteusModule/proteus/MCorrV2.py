from Transport import *
import cMCorrV2,cMCorr2D,cMCorrQ,cMCorr2DQ

class OneLevelMCorrV2(OneLevelTransport):
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
                 movingDomain=False):#,
        import pdb
#         OneLevelTransport.__init__(self,
#                                    uDict,
#                                    phiDict,
#                                    testSpaceDict,
#                                    matType,
#                                    dofBoundaryConditionsDict,
#                                    dofBoundaryConditionsSetterDict,
#                                    coefficients,
#                                    elementQuadrature,
#                                    elementBoundaryQuadrature,
#                                    fluxBoundaryConditionsDict,
#                                    advectiveFluxBoundaryConditionsSetterDict,
#                                    diffusiveFluxBoundaryConditionsSetterDictDict,
#                                    stabilization,
#                                    shockCapturing,
#                                    conservativeFluxDict,
#                                    numericalFluxType,
#                                    TimeIntegrationClass,
#                                    massLumping,
#                                    reactionLumping,
#                                    options,
#                                    name,
#                                    reuse_trial_and_test_quadrature,
#                                    sd,
#                                    movingDomain)
        import Comm
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
#         for ck,phi in phiDict.iteritems():
#             if coefficients.potential.has_key(ck):
#                 for cj in coefficients.potential[ck].keys():
#                     self.dphi[(ck,cj)] = FiniteElementFunction(phi.femSpace)
#             else:
#                 self.dphi[(ck,ck)] = FiniteElementFunction(phi.femSpace)
        #check for nonlinearities in the diffusion coefficient that don't match the potential
#         for ci,ckDict in coefficients.diffusion.iteritems():
#             #for ck,cjDict in coefficients.diffusion.iteritems(): #cek: bug?
#             for ck,cjDict in ckDict.iteritems():
#                 for cj in cjDict.keys():
#                     if not self.dphi.has_key((ck,cj)):
#                         self.dphi[(ck,cj)] = FiniteElementFunction(phi.femSpace)
        self.matType = matType
        #mwf try to reuse test and trial information across components if spaces are the same
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
        if type(self.u[0].femSpace) == C0_AffineLinearOnSimplexWithNodalBasis:
            print self.nQuadraturePoints_element
            if self.nSpace_global == 3:
                assert(self.nQuadraturePoints_element == 5)
            elif self.nSpace_global == 2:
                assert(self.nQuadraturePoints_element == 6)
            elif self.nSpace_global == 1:
                assert(self.nQuadraturePoints_element == 3)

            print self.nElementBoundaryQuadraturePoints_elementBoundary
            if self.nSpace_global == 3:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 2:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 1:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)

        #pdb.set_trace()
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
        self.q[('w*dV_r',0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0]),'d')
        self.q[('grad(w)*dV_a',0,0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0],self.nSpace_global),'d')
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('grad(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('r',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
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
        self.globalResidualDummy = None
    #mwf these are getting called by redistancing classes,
    def calculateCoefficients(self):
        pass
    def calculateElementResidual(self):
        if self.globalResidualDummy != None:
            self.getResidual(self.u[0].dof,self.globalResidualDummy)
    def getResidual(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.setUnknowns(u)
        #no flux boundary conditions
        self.elementResidual[0].fill(0.0)
        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
            cResidual = cMCorrV2.calculateResidual
            if self.nSpace_global == 2:
                cResidual = cMCorr2D.calculateResidual
            elif self.nSpace_global == 1:
                cResidual = cMCorr1D.calculateResidual
        elif isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            cResidual = cMCorrQ.calculateResidual
            if self.nSpace_global == 2:
                cResidual = cMCorr2DQ.calculateResidual
            elif self.nSpace_global == 1:
                cResidual = cMCorr1DQ.calculateResidual

        #cMCorr.calculateResidual(self.mesh.nElements_global,
        cResidual(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            #physics
            self.mesh.nElements_global,
                                 self.coefficients.epsFactHeaviside,
                                 self.coefficients.epsFactDirac,
                                 self.coefficients.epsFactDiffusion,
                                 self.u[0].femSpace.dofMap.l2g,
                                 self.mesh.elementDiametersArray,
                                 self.u[0].dof,
                                 self.q[('v',0)], 
                                 self.q[('grad(v)',0)], 
                                 self.q[('w*dV_r',0)], 
                                 self.q[('grad(w)*dV_a',0,0)], 
                                 self.coefficients.q_u_ls,
                                 self.coefficients.q_H_vof,
                                 self.q[('u',0)],
                                 self.q[('r',0)],
                                 self.elementResidual[0], 
                                 self.offset[0],self.stride[0],
                                 r)
        log("Global residual",level=9,data=r)
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy == None:
            self.globalResidualDummy = numpy.zeros(r.shape,'d')
    def getJacobian(self,jacobian):
        import superluWrappers
        import numpy
        import pdb
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,jacobian)

        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
            cJacobian = cMCorrV2.calculateJacobian
            if self.nSpace_global == 2:
                cJacobian = cMCorr2D.calculateJacobian
            elif self.nSpace_global == 1:
                cJacobian = cMCorr1D.calculateJacobian
        elif isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            cJacobian = cMCorrQ.calculateJacobian
            if self.nSpace_global == 2:
                cJacobian = cMCorr2DQ.calculateJacobian
            elif self.nSpace_global == 1:
                cJacobian = cMCorr1DQ.calculateJacobian
				      
        #cMCorr.calculateJacobian(self.mesh.nElements_global,
        cJacobian(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.nElements_global,
                                 self.coefficients.epsFactHeaviside,
                                 self.coefficients.epsFactDirac,
                                 self.coefficients.epsFactDiffusion,
                                 self.u[0].femSpace.dofMap.l2g,
                                 self.mesh.elementDiametersArray,
                                 self.u[0].dof,
                                 self.q[('v',0)], 
                                 self.q[('grad(v)',0)], 
                                 self.q[('w*dV_r',0)], 
                                 self.q[('grad(w)*dV_a',0,0)], 
                                 self.coefficients.q_u_ls,
                                 self.coefficients.q_H_vof,
                                 self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
                                 jacobian)
        log("Jacobian ",level=10,data=jacobian)
        #mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian
    
    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.
        
        This function should be called only when the mesh changes.
        """
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
	#
        #get physical locations of quadrature points and jacobian information there
	#assume all components live on the same mesh
        #
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
        cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('r',0)],
                                             self.q['abs(det(J))'],
                                             self.q[('w',0)],
                                             self.q[('w*dV_r',0)])
        self.testSpace[0].getBasisGradientValues(self.elementQuadraturePoints,
                                                  self.q['inverse(J)'],
                                                  self.q[('grad(w)',0)])
        cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('a',0,0)],
                                                      self.q['abs(det(J))'],
                                                      self.q[('grad(w)',0)],
                                                      self.q[('grad(w)*dV_a',0,0)])
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization != None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing != None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        pass
    def calculateExteriorElementBoundaryQuadrature(self):
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
    def estimate_mt(self):
        pass
