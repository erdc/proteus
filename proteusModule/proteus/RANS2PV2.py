from Transport import *
import cRANS2PV2,cRANS2P2D,cRANS2PQ#,cRANS2P2DQ

class OneLevelRANS2PV2(OneLevelTransport):
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
                 stressTraceBoundaryConditionsSetterDictDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='RANS2PV2',
                 reuse_trial_and_test_quadrature=True,
                 sd = True,
                 movingDomain=False):
        #
        #set the objects describing the method and boundary conditions
        #
        self.movingDomain=movingDomain
        self.tLast_mesh=None
        #
        #cek todo clean up these flags in the optimized version
        self.bcsTimeDependent=options.bcsTimeDependent
        self.bcsSet=False
        self.name=name
        self.sd=sd
        self.lowmem=True
        self.timeTerm=True#allow turning off  the  time derivative
        self.testIsTrial=True
        self.phiTrialIsTrial=True            
        self.u = uDict
        self.Hess=False
        if isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            self.Hess=True
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
        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
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
        #self.q['det(J)'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        #self.q['abs(det(J))'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        #self.q['J'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global,self.nSpace_global),'d')
        #self.q['inverse(J)'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global,self.nSpace_global),'d')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
	#self.ebqe['g'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
	#			       self.nElementBoundaryQuadraturePoints_elementBoundary,
        #max(1,self.nSpace_global-1),
        #max(1,self.nSpace_global-1)),
        #'d')
        #self.ebqe['inverse(J)'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global,self.nSpace_global),'d')
        #self.ebqe['hat(x)'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        #self.ebqe['bar(x)'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        #self.ebqe['sqrt(det(g))'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        #self.ebqe[('n')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #ebq for post-processing
        #self.ebq['x'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
	#self.ebq['g'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,
				      #  self.nElementBoundaryQuadraturePoints_elementBoundary,
				      #  max(1,self.nSpace_global-1),
				      #  max(1,self.nSpace_global-1)),
				      # 'd')
        #self.ebq['inverse(J)'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary,
                                              # self.nSpace_global,self.nSpace_global),'d')
        #self.ebq['hat(x)'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        #self.ebq['bar(x)'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        #self.ebq['sqrt(det(g))'] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        #self.ebq[('n')] = numpy.zeros((self.mesh.nElements_global,self.mesh.nElementBoundaries_element,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #ebq_global for post-processing
        #self.ebq_global['x'] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
	#self.ebq_global['g'] = numpy.zeros((self.mesh.nElementBoundaries_global,
        #                                    self.nElementBoundaryQuadraturePoints_elementBoundary,
        #                                    max(1,self.nSpace_global-1),
        #                                    max(1,self.nSpace_global-1)),
        #                                   'd')
        #self.ebq_global['inverse(J)'] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,
        #                                             self.nSpace_global,self.nSpace_global),'d')
        #self.ebq_global['hat(x)'] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        #self.ebq_global['bar(x)'] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        #self.ebq_global['sqrt(det(g))'] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        #self.ebq_global[('n')] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #need to calculate
        self.ebq_global[('totalFlux',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebq_global[('velocityAverage',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        #shape
        #self.q[('v',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0]),'d')
        #self.q[('grad(v)',0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[0],self.nSpace_global),'d')
        #self.q[('w*dV_r',0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1]),'d')
        #self.q[('grad(w)*dV_f',0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1],self.nSpace_global),'d')
        #self.ebqe[('v',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[0]),'d')
        #self.ebqe[('grad(v)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[0],self.nSpace_global),'d')
        #self.ebqe[('w*dS_f',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[0]),'d')
        #if self.nDOF_trial_element[1] != self.nDOF_trial_element[0]:
        #    self.q[('v',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1]),'d')
        #    self.q[('grad(v)',1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1],self.nSpace_global),'d')
        #    self.q[('w*dV_r',1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1]),'d')
        #    self.q[('grad(w)*dV_f',1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1],self.nSpace_global),'d')
        #    self.ebqe[('v',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[1]),'d')
        #    self.ebqe[('grad(v)',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[1],self.nSpace_global),'d')
        #self.ebqe[('w*dS_f',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nDOF_trial_element[1]),'d')
        #else:
        #    self.q[('v',1)] = self.q[('v',0)]
        #    self.q[('grad(v)',1)] = self.q[('grad(v)',0)]
        #    self.q[('w*dV_r',1)] =  self.q[('w*dV_r',0)]
        #    self.q[('grad(w)*dV_f',1)] = self.q[('grad(w)*dV_f',0)] 
        #    self.ebqe[('v',1)] = self.ebqe[('v',0)]
        #    self.ebqe[('grad(v)',1)] = self.ebqe[('grad(v)',0)]
        #    self.ebqe[('w*dS_f',1)] = self.ebqe[('w*dS_f',0)]
        #if self.Hess:
        #    self.q[('Hess(w)',1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1],self.nSpace_global,self.nSpace_global),'d')
        #    self.q[('Hess(w)*dV_a',1,1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nDOF_trial_element[1],self.nSpace_global,self.nSpace_global),'d')
        #else:
        #    #dummy
        #    self.q[('Hess(w)',1)] =  numpy.zeros((1,),'d')
        #    self.q[('Hess(w)*dV_a',1,1)] =  self.q[('Hess(w)',1)]
        #self.q[('Hess(v)',1)] =  self.q[('Hess(w)',1)]
        #self.q[('Hess(v)',2)] =  self.q[('Hess(w)',1)]
        #self.q[('Hess(v)',3)] =  self.q[('Hess(w)',1)]
        #self.q[('Hess(w)',2)] =  self.q[('Hess(w)',1)]
        #self.q[('Hess(w)',3)] =  self.q[('Hess(w)',1)]
        #self.q[('Hess(w)*dV_a',2,2)] =  self.q[('Hess(w)*dV_a',1,1)]
        #self.q[('Hess(w)*dV_a',3,3)] =  self.q[('Hess(w)*dV_a',1,1)]
        #self.q[('v',2)] = self.q[('v',1)]
        #self.q[('grad(v)',2)] = self.q[('grad(v)',1)]
        #self.q[('w*dV_r',2)] =  self.q[('w*dV_r',1)]
        #self.q[('grad(w)*dV_f',2)] = self.q[('grad(w)*dV_f',1)] 
        #self.ebqe[('v',2)] = self.ebqe[('v',1)]
        #self.ebqe[('grad(v)',2)] = self.ebqe[('grad(v)',1)]
        #self.ebqe[('w*dS_f',2)] = self.ebqe[('w*dS_f',1)]
        #self.q[('v',3)] = self.q[('v',1)]
        #self.q[('grad(v)',3)] = self.q[('grad(v)',1)]
        #self.q[('w*dV_r',3)] =  self.q[('w*dV_r',1)]
        #self.q[('grad(w)*dV_f',3)] = self.q[('grad(w)*dV_f',1)] 
        #self.ebqe[('v',3)] = self.ebqe[('v',1)]
        #self.ebqe[('grad(v)',3)] = self.ebqe[('grad(v)',1)]
        #self.ebqe[('w*dS_f',3)] = self.ebqe[('w*dS_f',1)]
        #for ci in range(self.nc):
        #    self.q[('w*dV_m',ci)] = self.q[('w*dV_r',ci)]
        #    self.q[('w',ci)] = self.q[('v',ci)]
        #    self.q[('grad(w)',ci)] = self.q[('grad(v)',ci)]
        #    self.ebqe[('w',ci)] = self.ebqe[('v',ci)]
        #    self.ebqe[('grad(w)',ci)] = self.ebqe[('grad(v)',ci)]
        #self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m',1)] = self.q[('u',1)]
        self.q[('m',2)] = self.q[('u',2)]
        self.q[('m',3)] = self.q[('u',3)]
        self.q[('m_last',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('f',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('velocity',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',1,1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',2,2)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',3,3)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',3,3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('penalty')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',3,3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('velocity',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
        #
        #show quadrature
        #
        log("Dumping quadrature shapes for model %s" % self.name,level=9)
        log("Element quadrature array (q)", level=9)
        for (k,v) in self.q.iteritems(): log(str((k,v.shape)),level=9)
        log("Element boundary quadrature (ebq)",level=9) 
        for (k,v) in self.ebq.iteritems(): log(str((k,v.shape)),level=9)
        log("Global element boundary quadrature (ebq_global)",level=9)
        for (k,v) in self.ebq_global.iteritems(): log(str((k,v.shape)),level=9)
        log("Exterior element boundary quadrature (ebqe)",level=9)
        for (k,v) in self.ebqe.iteritems(): log(str((k,v.shape)),level=9)
        log("Interpolation points for nonlinear diffusion potential (phi_ip)",level=9)
        for (k,v) in self.phi_ip.iteritems(): log(str((k,v.shape)),level=9)
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
        for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc',ci)].shape,'i')
            for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
            for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                self.ebqe[('diffusiveFlux_bc_flag',ck,ci)] = numpy.zeros(self.ebqe[('diffusiveFlux_bc',ck,ci)].shape,'i')
                for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                    self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        self.numericalFlux.setDirichletValues(self.ebqe)

        #cek/ido todo replace python loops in modules with optimized code if possible/necessary
        self.forceStrongConditions=True ##False
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(self.u[cj].femSpace,dofBoundaryConditionsSetterDict[cj],weakDirichletConditions=False)

    def getResidual(self,u,r):
        """
        Calculate the element residuals and add in to the global residual
        """
        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
        #if isinstance(self.u[0].femSpace,C0_AffineLinearOnCubeWithNodalBasis):
            cResidual = cRANS2PV2.calculateResidual
            if self.nSpace_global == 2:
                cResidual = cRANS2P2D.calculateResidual
        elif isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            cResidual = cRANS2PQ.calculateResidual
        #ido todo add cases for hex's and nurbs
            #if self.nSpace_global == 2:
            #    cResidual = cRANS2P2DQ.calculateResidual
        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        #cek todo put in logic to skip if BC's don't depend on t or u
        #hack
        if self.bcsTimeDependent or not self.bcsSet:
            self.bcsSet=True
            #Dirichlet boundary conditions
            self.numericalFlux.setDirichletValues(self.ebqe)
            #Flux boundary conditions
            for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
                for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                    if self.coefficients.advection.has_key(ci):
                        self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                        self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
                for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                    for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                        self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                        self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        r.fill(0.0)
        self.elementResidual[0].fill(0.0)
        self.elementResidual[1].fill(0.0) 
        self.elementResidual[2].fill(0.0) 
        self.elementResidual[3].fill(0.0)
        self.Ct_sge = 4.0
        self.Cd_sge = 144.0
 
        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                    self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN],self.timeIntegration.t)
 
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
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            #physics
            self.mesh.nElements_global,
            self.timeIntegration.alpha_bdf,
            self.coefficients.epsFact_density,
            self.coefficients.epsFact,
            self.coefficients.sigma,
            self.coefficients.rho_0,
            self.coefficients.nu_0,
            self.coefficients.rho_1,
            self.coefficients.nu_1,
            self.Ct_sge,
            self.Cd_sge,
            self.shockCapturing.shockCapturingFactor,
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.coefficients.g,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.timeIntegration.m_tmp[3],
            self.q[('f',0)],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.timeIntegration.beta_bdf[3],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.q[('numDiff',1,1)], 
            self.q[('numDiff',2,2)], 
            self.q[('numDiff',3,3)],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.shockCapturing.numDiff_last[3],
            self.elementResidual[0], 
            self.elementResidual[1], 
            self.elementResidual[2], 
            self.elementResidual[3],
            self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
            self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
            self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
            self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
            self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
            self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
            self.offset[0],self.offset[1],self.offset[2],self.offset[3],
            self.stride[0],self.stride[1],self.stride[2],self.stride[3],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_phi,
            self.coefficients.ebqe_n,
            self.coefficients.ebqe_kappa,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.numericalFlux.isDOFBoundary[3],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('advectiveFlux_bc_flag',3)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.ebqe[('diffusiveFlux_bc_flag',3,3)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.ebqe[('advectiveFlux_bc',3)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.numericalFlux.ebqe[('u',3)],
            self.ebqe[('diffusiveFlux_bc',3,3)],
            self.q[('velocity',0)],
            self.ebqe[('velocity',0)],
            self.ebq_global[('totalFlux',0)])

	if self.forceStrongConditions:#
	    for cj in range(len(self.dirichletConditionsForceDOF)):#
		for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                     r[self.offset[cj]+self.stride[cj]*dofN] = 0

        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        log("Global residual",level=9,data=r)
        #mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1
    def getJacobian(self,jacobian):
        import superluWrappers
        import numpy
        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
        #if isinstance(self.u[0].femSpace,C0_AffineLinearOnCubeWithNodalBasis):
            cJacobian = cRANS2PV2.calculateJacobian
            if self.nSpace_global == 2:
                cJacobian = cRANS2P2D.calculateJacobian
        elif isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            cJacobian = cRANS2PQ.calculateJacobian
            #if self.nSpace_global == 2:
            #    cJacobian = cRANS2P2DQ.calculateJacobian
        #ido todo add cases for hex/nurbs
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
				       jacobian)
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
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            #element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.nElements_global,
            self.timeIntegration.alpha_bdf,
            self.coefficients.epsFact_density,
            self.coefficients.epsFact,
            self.coefficients.sigma,
            self.coefficients.rho_0,
            self.coefficients.nu_0,
            self.coefficients.rho_1,
            self.coefficients.nu_1,
            self.Ct_sge,
            self.Cd_sge,
            self.shockCapturing.shockCapturingFactor,
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.coefficients.g,
            self.coefficients.q_phi,
            self.coefficients.q_n,
            self.coefficients.q_kappa,
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.timeIntegration.beta_bdf[3],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.shockCapturing.numDiff_last[3],
            self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
            self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
            self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
            self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
            self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
            self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            self.csrRowIndeces[(0,1)],self.csrColumnOffsets[(0,1)],
            self.csrRowIndeces[(0,2)],self.csrColumnOffsets[(0,2)],
            self.csrRowIndeces[(0,3)],self.csrColumnOffsets[(0,3)],
            self.csrRowIndeces[(1,0)],self.csrColumnOffsets[(1,0)],
            self.csrRowIndeces[(1,1)],self.csrColumnOffsets[(1,1)],
            self.csrRowIndeces[(1,2)],self.csrColumnOffsets[(1,2)],
            self.csrRowIndeces[(1,3)],self.csrColumnOffsets[(1,3)],
            self.csrRowIndeces[(2,0)],self.csrColumnOffsets[(2,0)],
            self.csrRowIndeces[(2,1)],self.csrColumnOffsets[(2,1)],
            self.csrRowIndeces[(2,2)],self.csrColumnOffsets[(2,2)],
            self.csrRowIndeces[(2,3)],self.csrColumnOffsets[(2,3)],
            self.csrRowIndeces[(3,0)],self.csrColumnOffsets[(3,0)],
            self.csrRowIndeces[(3,1)],self.csrColumnOffsets[(3,1)],
            self.csrRowIndeces[(3,2)],self.csrColumnOffsets[(3,2)],
            self.csrRowIndeces[(3,3)],self.csrColumnOffsets[(3,3)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_phi,
            self.coefficients.ebqe_n,
            self.coefficients.ebqe_kappa,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.numericalFlux.isDOFBoundary[3],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('advectiveFlux_bc_flag',3)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.ebqe[('diffusiveFlux_bc_flag',3,3)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.ebqe[('advectiveFlux_bc',3)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.numericalFlux.ebqe[('u',3)],
            self.ebqe[('diffusiveFlux_bc',3,3)],
            self.csrColumnOffsets_eb[(0,0)],
            self.csrColumnOffsets_eb[(0,1)],
            self.csrColumnOffsets_eb[(0,2)],
            self.csrColumnOffsets_eb[(0,3)],
            self.csrColumnOffsets_eb[(1,0)],
            self.csrColumnOffsets_eb[(1,1)],
            self.csrColumnOffsets_eb[(1,2)],
            self.csrColumnOffsets_eb[(1,3)],
            self.csrColumnOffsets_eb[(2,0)],
            self.csrColumnOffsets_eb[(2,1)],
            self.csrColumnOffsets_eb[(2,2)],
            self.csrColumnOffsets_eb[(2,3)],
            self.csrColumnOffsets_eb[(3,0)],
            self.csrColumnOffsets_eb[(3,1)],
            self.csrColumnOffsets_eb[(3,2)],
            self.csrColumnOffsets_eb[(3,3)])

        #Load the Dirichlet conditions directly into residual
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
        self.u[1].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
	#
        #get physical locations of quadrature points and jacobian information there
	#assume all components live on the same mesh
        #
        # self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
        #                                           self.q['x'])
        # if self.movingDomain: 
        #     if self.tLast_mesh != None:
        #         self.q['xt'][:]=self.q['x']
        #         self.q['xt']-=self.q['x_last']
        #         alpha = 1.0/(self.t_mesh - self.tLast_mesh)
        #         self.q['xt']*=alpha
        #     else:
        #         self.q['xt'][:]=0.0
        #     self.q['x_last'][:]=self.q['x']
        # self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
        #                                                  self.q['J'],
        #                                                  self.q['inverse(J)'],
        #                                                  self.q['det(J)'])
        # self.q['abs(det(J))']=numpy.absolute(self.q['det(J)'])
        # #
        # # get physical space integration weights
        # #
        # self.q['dV'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        # cfemIntegrals.calculateIntegrationWeights(self.q['abs(det(J))'],
        #                                           self.elementQuadratureWeights[('u',0)],
        #                                           self.q['dV'])
        # for ci in range(self.nc): self.q[('dV_u',ci)] = self.q['dV']
        # #
        # #get shape information at the quadrature points
        # #
        # self.testSpace[0].getBasisValues(self.elementQuadraturePoints,
        #                                  self.q[('w',0)])
        # cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('r',1)],
        #                                      self.q['abs(det(J))'],
        #                                      self.q[('w',1)],
        #                                      self.q[('w*dV_r',1)])
        # self.testSpace[1].getBasisGradientValues(self.elementQuadraturePoints,
        #                                           self.q['inverse(J)'],
        #                                           self.q[('grad(w)',1)])
        # cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('f',0)],
        #                                               self.q['abs(det(J))'],
        #                                               self.q[('grad(w)',0)],
        #                                               self.q[('grad(w)*dV_f',0)])
        # if self.Hess:
        #     self.testSpace[1].getBasisHessianValues(self.elementQuadraturePoints,
        #                                             self.q['inverse(J)'],
        #                                             self.q[('Hess(w)',1)])
        #     cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('a',1,1)],
        #                                                   self.q['abs(det(J))'],
        #                                                   self.q[('Hess(w)',1)],
        #                                                   self.q[('Hess(w)*dV_a',1,1)])
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization != None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing != None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on element boundaries.

        This function should be called only when the mesh changes.
        """
        #cek todo, a local calculation of the element boundary quadrature permuations
        #
        #get physical locations of element boundary quadrature points
        #
	#assume all components live on the same mesh
        # self.u[0].femSpace.elementMaps.getValuesTrace(self.elementBoundaryQuadraturePoints,
	# 					      self.ebq['x'])
        #
        #get metric tensor and unit normals
        #
        # if self.movingDomain:
        #     if self.tLast_mesh != None:
        #         self.ebq['xt'][:]=self.ebq['x']
        #         self.ebq['xt']-=self.ebq['x_last']
        #         alpha = 1.0/(self.t_mesh - self.tLast_mesh)
        #         self.ebq['xt']*=alpha
        #     else:
        #         self.ebq['xt'][:]=0.0
        #     self.ebq['x_last'][:]=self.ebq['x']
        #     self.u[0].femSpace.elementMaps.getJacobianValuesTrace_movingDomain(self.elementBoundaryQuadraturePoints,
        #                                                                        self.ebq['xt'],
        #                                                                        self.ebq['inverse(J)'],
        #                                                                        self.ebq['g'],
        #                                                                        self.ebq['sqrt(det(g))'],
        #                                                                        self.ebq['n'])
        # else:
        #     self.u[0].femSpace.elementMaps.getJacobianValuesTrace(self.elementBoundaryQuadraturePoints,
        #                                                           self.ebq['inverse(J)'],
        #                                                           self.ebq['g'],
        #                                                           self.ebq['sqrt(det(g))'],
        #                                                           self.ebq['n'])

        # cfemIntegrals.copyLeftElementBoundaryInfo(self.mesh.elementBoundaryElementsArray,
        #                                           self.mesh.elementBoundaryLocalElementBoundariesArray,
        #                                           self.mesh.exteriorElementBoundariesArray,
        #                                           self.mesh.interiorElementBoundariesArray,
        #                                           self.ebq['x'],
        #                                           self.ebq['n'],
        #                                           self.ebq_global['x'],
        #                                           self.ebq_global['n'])
        # if self.movingDomain:
        #     cfemIntegrals.copyLeftElementBoundaryInfo_movingDomain(self.mesh.elementBoundaryElementsArray,
        #                                                            self.mesh.elementBoundaryLocalElementBoundariesArray,
        #                                                            self.mesh.exteriorElementBoundariesArray,
        #                                                            self.mesh.interiorElementBoundariesArray,
        #                                                            self.ebq['xt'])
        # #now map the physical points back to the reference element
        # #assume all components live  on same mesh
        # self.u[0].femSpace.elementMaps.getInverseValuesTrace(self.ebq['inverse(J)'],self.ebq['x'],self.ebq['hat(x)'])
        # self.u[0].femSpace.elementMaps.getPermutations(self.ebq['hat(x)'])
        #
        #since the points on the reference boundary may be reordered on many right element boundaries, we
        #have to use an array of reference boundary points on all element boundaries
        #first copy the left reference element boundary quadrature points from the reference element boundary
        #
        #get the shape information at the reference element boundary quadrature points 
	#
        # self.testSpace[0].getBasisValuesTrace(self.u[0].femSpace.elementMaps.permutations,
        #                                        self.ebq['hat(x)'],
        #                                        self.ebq[('w',0)])
        # cfemIntegrals.calculateWeightedShapeTrace(self.elementBoundaryQuadratureWeights[('u',0)],
        #                                           self.ebq['sqrt(det(g))'],
        #                                           self.ebq[('w',0)],
        #                                           self.ebq[('w*dS_u',0)])
        # self.u[0].femSpace.getBasisGradientValuesTrace(self.u[0].femSpace.elementMaps.permutations,
        #                                                 self.ebq['hat(x)'],
        #                                                 self.ebq['inverse(J)'],
        #                                                 self.ebq[('grad(v)',0)])
        # cfemIntegrals.calculateElementBoundaryIntegrationWeights(self.ebq['sqrt(det(g))'],
        #                                                         self.elementBoundaryQuadratureWeights[('u',0)],
        #                                                         self.ebq[('dS_u',0)])
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
        self.u[1].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        # #
        # #get metric tensor and unit normals
        # #
        # if self.movingDomain:
        #     if self.tLast_mesh != None:
        #         self.ebqe['xt'][:]=self.ebqe['x']
        #         self.ebqe['xt']-=self.ebqe['x_last']
        #         alpha = 1.0/(self.t_mesh - self.tLast_mesh)
        #         self.ebqe['xt']*=alpha
        #     else:
        #         self.ebqe['xt'][:]=0.0
        #     self.ebqe['x_last'][:]=self.ebqe['x']
        #     self.u[0].femSpace.elementMaps.getJacobianValuesGlobalExteriorTrace_movingDomain(self.elementBoundaryQuadraturePoints,
        #                                                                                      self.ebqe['xt'],
        #                                                                                      self.ebqe['inverse(J)'],
        #                                                                                      self.ebqe['g'],
        #                                                                                      self.ebqe['sqrt(det(g))'],
        #                                                                                      self.ebqe['n'])
        # else:
        #     self.u[0].femSpace.elementMaps.getJacobianValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
        #                                                                         self.ebqe['inverse(J)'],
        #                                                                         self.ebqe['g'],
        #                                                                         self.ebqe['sqrt(det(g))'],
        #                                                                         self.ebqe['n'])
        # #now map the physical points back to the reference element
        # #assume all components live  on same mesh
        # self.u[0].femSpace.elementMaps.getInverseValuesGlobalExteriorTrace(self.ebqe['inverse(J)'],self.ebqe['x'],self.ebqe['hat(x)'])
        # #
        # #since the points on the reference boundary may be reordered on many right element boundaries, we
        # #have to use an array of reference boundary points on all element boundaries
        # #first copy the left reference element boundary quadrature points from the reference element boundary
        # self.testSpace[0].getBasisValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
        #                                                      self.ebqe[('w',0)])
        # cfemIntegrals.calculateWeightedShapeGlobalExteriorTrace(self.mesh.exteriorElementBoundariesArray,
        #                                                         self.mesh.elementBoundaryElementsArray,
        #                                                         self.mesh.elementBoundaryLocalElementBoundariesArray,
        #                                                         self.elementBoundaryQuadratureWeights[('f',0)],
        #                                                         self.ebqe['sqrt(det(g))'],
        #                                                         self.ebqe[('w',0)],
        #                                                         self.ebqe[('w*dS_f',0)])
        # self.u[0].femSpace.getBasisGradientValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
        #                                                               self.ebqe['inverse(J)'],
        #                                                               self.ebqe[('grad(v)',0)])
        # #setup flux boundary conditions
        self.fluxBoundaryConditionsObjectsDict = dict([(cj,FluxBoundaryConditions(self.mesh,
                                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                  self.ebqe[('x')],
                                                                                  self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        # self.ebqe['dS'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        # cfemIntegrals.calculateIntegrationWeights(self.ebqe['sqrt(det(g))'],
        #                                           self.elementBoundaryQuadratureWeights[('u',0)],
        #                                           self.ebqe['dS'])
        # for ci in range(self.nc): self.ebqe[('dS_u',ci)] = self.ebqe['dS']
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
    def estimate_mt(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
            cva = cRANS2PV2.calculateVelocityAverage
            if self.nSpace_global == 2:
                cva = cRANS2P2D.calculateVelocityAverage
        elif isinstance(self.u[0].femSpace,C0_AffineQuadraticOnSimplexWithNodalBasis):
            cva = cRANS2PQ.calculateVelocityAverage
            #if self.nSpace_global == 2:
            #    cva = cRANS2P2DQ.calculateVelocityAverage
        #cek todo memory optimized version
        # cva(self.u[1].femSpace.elementMaps.permutations,
        #     self.mesh.nExteriorElementBoundaries_global,
        #     self.mesh.exteriorElementBoundariesArray,
        #     self.mesh.nInteriorElementBoundaries_global,
        #     self.mesh.interiorElementBoundariesArray,
        #     self.mesh.elementBoundaryElementsArray,
        #     self.mesh.elementBoundaryLocalElementBoundariesArray,
        #     self.u[1].femSpace.dofMap.l2g,
        #     self.u[1].dof,
        #     self.u[2].dof,
        #     self.u[3].dof,
        #     self.ebq[('v',0)], 
        #     self.ebqe[('velocity',0)],
        #     self.ebq_global[('velocityAverage',0)])
        OneLevelTransport.calculateAuxiliaryQuantitiesAfterStep(self)
