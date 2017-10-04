import proteus
from .cRichards import *

class Coefficients(proteus.TransportCoefficients.TC_base):
    """
    version of Re where element material type id's used in evals
    """
    from proteus.ctransportCoefficients import conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2
    from proteus.ctransportCoefficients import conservativeHeadRichardsMualemVanGenuchten_sd_het
    def __init__(self,
                 nd,
                 Ksw_types,
                 vgm_n_types,
                 vgm_alpha_types,
                 thetaR_types,
                 thetaSR_types,
                 gravity,
                 density,
                 beta,
                 diagonal_conductivity=True,
                 getSeepageFace=None):
        variableNames=['pressure_head']
        nc=1
        mass={0:{0:'nonlinear'}}
        advection={0:{0:'nonlinear'}}
        diffusion={0:{0:{0:'nonlinear'}}}
        potential={0:{0:'u'}}
        reaction={0:{0:'linear'}}
        hamiltonian={}
        self.getSeepageFace=getSeepageFace
        self.gravity=gravity
        self.rho = density
        self.beta=beta
        self.vgm_n_types = vgm_n_types
        self.vgm_alpha_types = vgm_alpha_types
        self.thetaR_types    = thetaR_types
        self.thetaSR_types   = thetaSR_types
        self.elementMaterialTypes = None
        self.exteriorElementBoundaryTypes  = None
        self.materialTypes_q    = None
        self.materialTypes_ebq  = None
        self.materialTypes_ebqe  = None
        self.nd = nd
        self.nMaterialTypes = len(thetaR_types)
        self.q = {}; self.ebqe = {}; self.ebq = {}; self.ebq_global={}
        #try to allow some flexibility in input of permeability/conductivity tensor
        self.diagonal_conductivity = diagonal_conductivity
        self.Ksw_types_in = Ksw_types
        if self.diagonal_conductivity:
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd+1,dtype='i'),
                                             numpy.arange(self.nd,dtype='i'))}

            assert len(Ksw_types.shape) in [1,2], "if diagonal conductivity true then Ksw_types scalar or vector of diagonal entries"
            #allow scalar input Ks
            if len(Ksw_types.shape)==1:
                self.Ksw_types = numpy.zeros((self.nMaterialTypes,self.nd),'d')
                for I in range(self.nd):
                    self.Ksw_types[:,I] = Ksw_types
            else:
                self.Ksw_types = Ksw_types
        else: #full
            sparseDiffusionTensors = {(0,0):(numpy.arange(self.nd**2+1,step=self.nd,dtype='i'),
                                             numpy.array([range(self.nd) for row in range(self.nd)],dtype='i'))}
            assert len(Ksw_types.shape) in [1,2], "if full tensor conductivity true then Ksw_types scalar or 'flattened' row-major representation of entries"
            if len(Ksw_types.shape)==1:
                self.Ksw_types = numpy.zeros((self.nMaterialTypes,self.nd**2),'d')
                for I in range(self.nd):
                    self.Ksw_types[:,I*self.nd+I] = Ksw_types
            else:
                assert Ksw_types.shape[1] == self.nd**2
                self.Ksw_types = Ksw_types
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         variableNames,
                         sparseDiffusionTensors = sparseDiffusionTensors,
                         useSparseDiffusion = True)
            
    def initializeMesh(self,mesh):
        from proteus.SubsurfaceTransportCoefficients import BlockHeterogeneousCoefficients
        #import pdb
        #pdb.set_trace()
        self.elementMaterialTypes,self.exteriorElementBoundaryTypes,self.elementBoundaryTypes = BlockHeterogeneousCoefficients(mesh).initializeMaterialTypes()
        #want element boundary material types for evaluating heterogeneity
        #not boundary conditions
        self.isSeepageFace = numpy.zeros((mesh.nExteriorElementBoundaries_global),'i')
        if self.getSeepageFace != None:
            for ebNE in range(mesh.nExteriorElementBoundaries_global):
                #mwf missing ebNE-->ebN?
                ebN = mesh.exteriorElementBoundariesArray[ebNE]
                #print "eb flag",mesh.elementBoundaryMaterialTypes[ebN]
                #print self.getSeepageFace(mesh.elementBoundaryMaterialTypes[ebN])
                self.isSeepageFace[ebNE] = self.getSeepageFace(mesh.elementBoundaryMaterialTypes[ebN])
        #print self.isSeepageFace
    def initializeElementQuadrature(self,t,cq):
        self.materialTypes_q = self.elementMaterialTypes
        self.q_shape = cq[('u',0)].shape
#        cq['Ks'] = numpy.zeros(self.q_shape,'d')
#        for k in range(self.q_shape[1]):
#            cq['Ks'][:,k] = self.Ksw_types[self.elementMaterialTypes,0]
        self.q[('vol_frac',0)] = numpy.zeros(self.q_shape,'d')
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.materialTypes_ebq = numpy.zeros(cebq[('u',0)].shape[0:2],'i')
        self.ebq_shape = cebq[('u',0)].shape
        for ebN_local in range(self.ebq_shape[1]):
            self.materialTypes_ebq[:,ebN_local] = self.elementMaterialTypes
        self.ebq[('vol_frac',0)] = numpy.zeros(self.ebq_shape,'d')

    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        self.materialTypes_ebqe = self.exteriorElementBoundaryTypes
        self.ebqe_shape = cebqe[('u',0)].shape
        self.ebqe[('vol_frac',0)] = numpy.zeros(self.ebqe_shape,'d')
        #
    def evaluate(self,t,c):
        if c[('u',0)].shape == self.q_shape:
            materialTypes = self.materialTypes_q
            vol_frac = self.q[('vol_frac',0)]
        elif c[('u',0)].shape == self.ebqe_shape:
            materialTypes = self.materialTypes_ebqe
            vol_frac = self.ebqe[('vol_frac',0)]
        elif c[('u',0)].shape == self.ebq_shape:
            materialTypes = self.materialTypes_ebq
            vol_frac = self.ebq[('vol_frac',0)]
        else:
            assert False, "no materialType found to match c[('u',0)].shape= %s " % c[('u',0)].shape
        self.conservativeHeadRichardsMualemVanGenuchten_sd_het(self.sdInfo[(0,0)][0],
                                                               self.sdInfo[(0,0)][1],
                                                               materialTypes,
                                                               self.rho,
                                                               self.beta,
                                                               self.gravity,
                                                               self.vgm_alpha_types,
                                                               self.vgm_n_types,
                                                               self.thetaR_types,
                                                               self.thetaSR_types,
                                                               self.Ksw_types,
                                                               c[('u',0)],
                                                               c[('m',0)],
                                                               c[('dm',0,0)],
                                                               c[('f',0)],
                                                               c[('df',0,0)],
                                                               c[('a',0,0)],
                                                               c[('da',0,0,0)],
                                                               vol_frac)
        # print "Picard---------------------------------------------------------------"
        # c[('df',0,0)][:] = 0.0
        # c[('da',0,0,0)][:] = 0.0
#         self.conservativeHeadRichardsMualemVanGenuchtenHetEvaluateV2(materialTypes,
#                                                                      self.rho,
#                                                                      self.beta,
#                                                                      self.gravity,
#                                                                      self.vgm_alpha_types,
#                                                                      self.vgm_n_types,
#                                                                      self.thetaR_types,
#                                                                      self.thetaSR_types,
#                                                                      self.Ksw_types,
#                                                                      c[('u',0)],
#                                                                      c[('m',0)],
#                                                                      c[('dm',0,0)],
#                                                                      c[('f',0)],
#                                                                      c[('df',0,0)],
#                                                                      c[('a',0,0)],
#                                                                      c[('da',0,0,0)])
        #mwf debug
        if (numpy.isnan(c[('da',0,0,0)]).any() or
            numpy.isnan(c[('a',0,0)]).any() or
            numpy.isnan(c[('df',0,0)]).any() or
            numpy.isnan(c[('f',0)]).any() or
            numpy.isnan(c[('u',0)]).any() or
            numpy.isnan(c[('m',0)]).any() or
            numpy.isnan(c[('dm',0,0)]).any()):
            import pdb
            pdb.set_trace()

#         #mwf debug
#         if c[('u',0)].shape == self.q_shape:
#             c[('visPerm',0)]=c[('a',0,0)][:,:,0,0]

class LevelModel(proteus.Transport.OneLevelTransport):
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
                 movingDomain=False,
                 bdyNullSpace=False):
        self.bdyNullSpace=bdyNullSpace
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
#        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
#            print self.nQuadraturePoints_element
#            if self.nSpace_global == 3:
#                assert(self.nQuadraturePoints_element == 5)
#            elif self.nSpace_global == 2:
#                assert(self.nQuadraturePoints_element == 6)
#            elif self.nSpace_global == 1:
#                assert(self.nQuadraturePoints_element == 3)
#
#            print self.nElementBoundaryQuadraturePoints_elementBoundary
#            if self.nSpace_global == 3:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 2:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 1:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)

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
        #self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('grad(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q['velocity'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('m',0)] = self.q[('u',0)].copy()
        self.q[('mt',0)] = self.q[('u',0)].copy()
        self.q[('m_last',0)] = self.q[('u',0)].copy()
        self.q[('m_tmp',0)] = self.q[('u',0)].copy()
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',0,0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('grad(u)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe['velocity'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('penalty')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
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
        logEvent("Updating local to global mappings",2)
        self.updateLocal2Global()
        logEvent("Building time integration object",2)
        logEvent(memory("inflowBC, internalNodes,updateLocal2Global","OneLevelTransport"),level=4)
        #mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(self,integrateInterpolationPoints=True)
        else:
             self.timeIntegration = TimeIntegrationClass(self)
           
        if options != None:
            self.timeIntegration.setFromOptions(options)
        logEvent(memory("TimeIntegration","OneLevelTransport"),level=4)
        logEvent("Calculating numerical quadrature formulas",2)
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
        logEvent(memory("stride+offset","OneLevelTransport"),level=4)
        if numericalFluxType != None:
            if options is None or options.periodicDirichletConditions is None:
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
        logEvent(memory("numericalFlux","OneLevelTransport"),level=4)
        self.elementEffectiveDiametersArray  = self.mesh.elementInnerDiametersArray
        #use post processing tools to get conservative fluxes, None by default
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)  
        logEvent(memory("velocity postprocessor","OneLevelTransport"),level=4)
        #helper for writing out data storage
        from proteus import Archiver
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
            self.numericalFlux.isDOFBoundary = {0:numpy.zeros(self.ebqe[('u',0)].shape,'i')}
        if not hasattr(self.numericalFlux,'ebqe'):
            self.numericalFlux.ebqe = {('u',0):numpy.zeros(self.ebqe[('u',0)].shape,'d')}
        #TODO how to handle redistancing calls for calculateCoefficients,calculateElementResidual etc
        self.globalResidualDummy = None
        compKernelFlag=0
        self.vof = cRichards_base(self.nSpace_global,
                             self.nQuadraturePoints_element,
                             self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                             self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                             self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                             compKernelFlag)
        if self.movingDomain:
            self.MOVING_DOMAIN=1.0
        else:
            self.MOVING_DOMAIN=0.0
        #cek hack
        self.movingDomain=False
        self.MOVING_DOMAIN=0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape,'d')        
        self.forceStrongConditions=False
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(self.u[cj].femSpace,dofBoundaryConditionsSetterDict[cj],weakDirichletConditions=False)
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
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        #cek can put in logic to skip of BC's don't depend on t or u
        #Dirichlet boundary conditions
        self.numericalFlux.setDirichletValues(self.ebqe)
        #flux boundary conditions
        #cek hack, just using advective flux for flux BC for now
        for t,g in self.fluxBoundaryConditionsObjectsDict[0].advectiveFluxBoundaryConditionsDict.iteritems():
            self.ebqe[('advectiveFlux_bc',0)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
            self.ebqe[('advectiveFlux_bc_flag',0)][t[0],t[1]] = 1
        # for t,g in self.fluxBoundaryConditionsObjectsDict[0].diffusiveFluxBoundaryConditionsDict.iteritems():
        #     self.ebqe[('diffusiveFlux_bc',0)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
        #     self.ebqe[('diffusiveFlux_bc_flag',0)][t[0],t[1]] = 1
        #self.shockCapturing.lag=True
        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                    self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN],self.timeIntegration.t)
        self.vof.calculateResidual(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
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
            self.ebqe['penalty'],#double* ebqe_penalty_ext,
            self.mesh.elementMaterialTypes,#int* elementMaterialTypes,	
            self.coefficients.isSeepageFace,
            self.coefficients.sdInfo[(0,0)][0],#int* a_rowptr,
            self.coefficients.sdInfo[(0,0)][1],#int* a_colind,
            self.coefficients.rho,#double rho,
            self.coefficients.beta,#double beta,
            self.coefficients.gravity,#double* gravity,
            self.coefficients.vgm_alpha_types,#double* alpha,
            self.coefficients.vgm_n_types,#double* n,
            self.coefficients.thetaR_types,#double* thetaR,
            self.coefficients.thetaSR_types,#double* thetaSR,
            self.coefficients.Ksw_types,#double* KWs,            
	    False,#self.coefficients.useMetrics, 
            self.timeIntegration.alpha_bdf,
            0,#self.shockCapturing.lag,
            0.0,#cek hack self.shockCapturing.shockCapturingFactor,
	    0.0,#self.coefficients.sc_uref, 
	    0.0,#self.coefficients.sc_beta,
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.u[0].dof,#cek hack
            self.q['velocity'],#self.coefficients.q_v,
            self.timeIntegration.m_tmp[0],
            self.q[('u',0)],
            self.timeIntegration.beta_bdf[0],
            self.q[('cfl',0)],
            self.q[('cfl',0)],#cek hack self.shockCapturing.numDiff[0],
            self.q[('cfl',0)],#cek hack self.shockCapturing.numDiff_last[0],
            self.offset[0],self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.ebqe['velocity'],#self.coefficients.ebqe_v,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('u',0)],#cek hack            self.coefficients.ebqe_phi,
            0.0,#cek hack self.coefficients.epsFact,
            self.ebqe[('u',0)],
            self.ebqe[('advectiveFlux',0)])
        #self.q[('mt',0)][:] =self.timeIntegration.m_tmp[0]
        #self.q[('mt',0)] *= self.timeIntegration.alpha_bdf
        #self.q[('mt',0)] += self.timeIntegration.beta_bdf[0]
        #self.timeIntegration.calculateElementCoefficients(self.q)
	if self.forceStrongConditions:#
	    for cj in range(len(self.dirichletConditionsForceDOF)):#
		for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                     r[self.offset[cj]+self.stride[cj]*dofN] = 0
        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        logEvent("Global residual",level=9,data=r)
        #mwf debug
        #pdb.set_trace()
        #mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy is None:
            self.globalResidualDummy = numpy.zeros(r.shape,'d')
    def getJacobian(self,jacobian):
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
				       jacobian)
        self.vof.calculateJacobian(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
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
            self.ebqe['penalty'],#double* ebqe_penalty_ext,
            self.mesh.elementMaterialTypes,#int* elementMaterialTypes,	
            self.coefficients.isSeepageFace,
            self.coefficients.sdInfo[(0,0)][0],#int* a_rowptr,
            self.coefficients.sdInfo[(0,0)][1],#int* a_colind,
            self.coefficients.rho,#double rho,
            self.coefficients.beta,#double beta,
            self.coefficients.gravity,#double* gravity,
            self.coefficients.vgm_alpha_types,#double* alpha,
            self.coefficients.vgm_n_types,#double* n,
            self.coefficients.thetaR_types,#double* thetaR,
            self.coefficients.thetaSR_types,#double* thetaSR,
            self.coefficients.Ksw_types,#double* KWs,            
	    False,#self.coefficients.useMetrics, 
            self.timeIntegration.alpha_bdf,
            0,#self.shockCapturing.lag,
            0.0,#self.shockCapturing.shockCapturingFactor,
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.q['velocity'],#self.coefficients.q_v,
            self.timeIntegration.beta_bdf[0],
            self.q[('cfl',0)],
            self.q[('cfl',0)],#cek hack self.shockCapturing.numDiff_last[0],
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.ebqe['velocity'],#self.coefficients.ebqe_v,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.csrColumnOffsets_eb[(0,0)])
        if self.forceStrongConditions:
            scaling = 1.0#probably want to add some scaling to match non-dirichlet diagonals in linear system 
            for cj in range(self.nc):
                for dofN in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.keys():
                    global_dofN = self.offset[cj]+self.stride[cj]*dofN
                    for i in range(self.rowptr[global_dofN],self.rowptr[global_dofN+1]):
                        if (self.colind[i] == global_dofN):
                            #print "RBLES forcing residual cj = %s dofN= %s global_dofN= %s was self.nzval[i]= %s now =%s " % (cj,dofN,global_dofN,self.nzval[i],scaling)
                            self.nzval[i] = scaling
                        else:
                            self.nzval[i] = 0.0
                            #print "RBLES zeroing residual cj = %s dofN= %s global_dofN= %s " % (cj,dofN,global_dofN)
        logEvent("Jacobian ",level=10,data=jacobian)
        #mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian
    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.
        
        This function should be called only when the mesh changes.
        """
        #self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
        #                                         self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization != None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing != None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        pass
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
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        self.fluxBoundaryConditionsObjectsDict = dict([(cj,FluxBoundaryConditions(self.mesh,
                                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                  self.ebqe[('x')],
                                                                                  getAdvectiveFluxBoundaryConditions=self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  getDiffusiveFluxBoundaryConditions=self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
    def estimate_mt(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass
