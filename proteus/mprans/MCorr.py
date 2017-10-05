import proteus
from proteus.mprans.cMCorr import *

class Coefficients(proteus.TransportCoefficients.TC_base):
    from proteus.ctransportCoefficients import levelSetConservationCoefficientsEvaluate
    from proteus.ctransportCoefficients import levelSetConservationCoefficientsEvaluate_sd
    def __init__(self,
                 applyCorrection=True,
                 epsFactHeaviside=0.0,
                 epsFactDirac=1.0,
                 epsFactDiffusion=2.0,
                 LSModel_index=3,
                 V_model=2,
                 me_model=5,
                 VOFModel_index=4,
                 checkMass=True,
                 sd=True,
                 nd=None,
                 applyCorrectionToDOF=True,
                 useMetrics=0.0,
                 useConstantH=False,
                 # mql. For edge based stabilization methods
                 useQuadraticRegularization=False,
                 edgeBasedStabilizationMethods=False):

        self.useQuadraticRegularization=useQuadraticRegularization
        self.edgeBasedStabilizationMethods=edgeBasedStabilizationMethods                 
        self.useConstantH=useConstantH
        self.useMetrics=useMetrics
        self.sd=sd
        self.checkMass=checkMass
        self.variableNames=['phiCorr']
        nc=1
        mass={}
        advection={}
        hamiltonian={}
        diffusion={0:{0:{0:'constant'}}}
        potential={0:{0:'u'}}
        reaction={0:{0:'nonlinear'}}
        #reaction={}
        if self.sd:
            assert nd is not None,"You must set the number of dimensions to use sparse diffusion in LevelSetConservationCoefficients"
            sdInfo = {(0,0):(numpy.arange(start=0,stop=nd+1,step=1,dtype='i'),
                             numpy.arange(start=0,stop=nd,step=1,dtype='i'))}
        else:
            sdInfo={}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         self.variableNames,
                         sparseDiffusionTensors=sdInfo,
                         useSparseDiffusion = sd)
        self.useQuadraticRegularization=useQuadraticRegularization
        self.edgeBasedStabilizationMethods=edgeBasedStabilizationMethods
        self.levelSetModelIndex=LSModel_index
        self.flowModelIndex=V_model
        self.epsFactHeaviside=epsFactHeaviside
        self.epsFactDirac=epsFactDirac
        self.epsFactDiffusion=epsFactDiffusion
        self.me_model=me_model
        self.VOFModelIndex=VOFModel_index
        self.useC = True
        self.applyCorrection=applyCorrection
        if self.applyCorrection:
            self.applyCorrectionToDOF = applyCorrectionToDOF
        else:
            self.applyCorrectionToDOF = False
        self.massConservationError=0.0
    def initializeMesh(self,mesh):
        self.h=mesh.h
        self.epsHeaviside = self.epsFactHeaviside*mesh.h
        self.epsDirac = self.epsFactDirac*mesh.h
        self.epsDiffusion = (self.epsFactDiffusion*mesh.h*
                             (mesh.h if self.useQuadraticRegularization==True else 1.))
    def attachModels(self,modelList):
        import copy
        logEvent("Attaching models in LevelSetConservation")
        #level set
        self.lsModel = modelList[self.levelSetModelIndex]
        self.q_u_ls  = modelList[self.levelSetModelIndex].q[('u',0)]
	self.q_n_ls  = modelList[self.levelSetModelIndex].q[('grad(u)',0)]

        self.ebqe_u_ls = modelList[self.levelSetModelIndex].ebqe[('u',0)]
        self.ebqe_n_ls = modelList[self.levelSetModelIndex].ebqe[('grad(u)',0)]

        if modelList[self.levelSetModelIndex].ebq.has_key(('u',0)):
            self.ebq_u_ls = modelList[self.levelSetModelIndex].ebq[('u',0)]
        else:
            self.ebq_u_ls = None

        #volume of fluid
        self.vofModel = modelList[self.VOFModelIndex]
        self.q_H_vof = modelList[self.VOFModelIndex].q[('u',0)]
        self.q_porosity = modelList[self.VOFModelIndex].coefficients.q_porosity
        self.ebqe_H_vof = modelList[self.VOFModelIndex].ebqe[('u',0)]
        if modelList[self.VOFModelIndex].ebq.has_key(('u',0)):
            self.ebq_H_vof = modelList[self.VOFModelIndex].ebq[('u',0)]
        else:
            self.ebq_H_vof = None
        #correction
        self.massCorrModel = modelList[self.me_model]
        self.massCorrModel.setMassQuadrature()
        self.vofModel.q[('m_last',0)][:] = self.vofModel.q[('m',0)]
        if self.checkMass:
            self.m_tmp = copy.deepcopy(self.massCorrModel.q[('r',0)])
            if self.checkMass:
                # self.vofGlobalMass = Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                #                                                 self.vofModel.q[('u',0)],
                #                                                 self.massCorrModel.mesh.nElements_owned)
                # self.lsGlobalMass = Norms.scalarHeavisideDomainIntegral(self.vofModel.q['dV'],
                #                                                         self.lsModel.q[('u',0)],
                #                                                         self.massCorrModel.mesh.nElements_owned)
                #self.vofGlobalMass = 0.0
                #self.lsGlobalMass = self.massCorrModel.calculateMass(self.lsModel.q[('u',0)])
                #logEvent("Attach Models MCorr: mass correction %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                #                                                                                self.massCorrModel.q[('r',0)],
                #                                                                                self.massCorrModel.mesh.nElements_owned),),level=2)
                self.fluxGlobal = 0.0
                self.totalFluxGlobal = 0.0
                self.vofGlobalMassArray = []#self.vofGlobalMass]
                self.lsGlobalMassArray = []#self.lsGlobalMass]
                self.vofGlobalMassErrorArray = []#self.vofGlobalMass - self.vofGlobalMassArray[0]]# + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral]
                self.lsGlobalMassErrorArray = []#self.lsGlobalMass - self.lsGlobalMassArray[0]]# + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral]
                self.fluxArray = []#0.0]#self.vofModel.coefficients.fluxIntegral]
                self.timeArray = []#self.vofModel.timeIntegration.t]
                #logEvent("Attach Models MCorr: Phase 0 mass after mass correction (VOF) %12.5e" % (self.vofGlobalMass,),level=2)
                #logEvent("Attach Models MCorr: Phase 0 mass after mass correction (LS) %12.5e" % (self.lsGlobalMass,),level=2)
                #logEvent("Attach Models MCorr: Phase  0 mass conservation (VOF) after step = %12.5e" % (self.vofGlobalMass - self.vofModel.coefficients.m_pre + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral,),level=2)
                #logEvent("Attach Models MCorr: Phase  0 mass conservation (LS) after step = %12.5e" % (self.lsGlobalMass - self.lsModel.coefficients.m_pre + self.vofModel.timeIntegration.dt*self.vofModel.coefficients.fluxIntegral,),level=2)
    def initializeElementQuadrature(self,t,cq):
        if self.sd and cq.has_key(('a',0,0)):
            cq[('a',0,0)].fill(self.epsDiffusion)
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        if self.sd and cebq.has_key(('a',0,0)):
            cebq[('a',0,0)].fill(self.epsDiffusion)
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        if self.sd and cebqe.has_key(('a',0,0)):
            cebqe[('a',0,0)].fill(self.epsDiffusion)
    def preStep(self,t,firstStep=False):
        if self.checkMass:
            logEvent("Phase 0 mass before mass correction (VOF) %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                                 self.vofModel.q[('m',0)],
                                                                                                 self.massCorrModel.mesh.nElements_owned),),level=2)
            logEvent("Phase 0 mass (primitive) before mass correction (LS) %12.5e" % (Norms.scalarSmoothedHeavisideDomainIntegral(self.epsFactHeaviside,
                                                                                                                 self.massCorrModel.elementDiameter,
                                                                                                                 self.vofModel.q['dV'],
                                                                                                                 self.lsModel.q[('m',0)],
                                                                                                                 self.massCorrModel.mesh.nElements_owned),),level=2)
            logEvent("Phase 0 mass (consistent) before mass correction (LS) %12.5e" % (self.massCorrModel.calculateMass(self.lsModel.q[('m',0)]),),level=2)
        copyInstructions = {'clear_uList':True}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        if self.applyCorrection:
            #ls
            self.lsModel.u[0].dof += self.massCorrModel.u[0].dof
            self.lsModel.q[('u',0)] += self.massCorrModel.q[('u',0)]
            self.lsModel.ebqe[('u',0)] += self.massCorrModel.ebqe[('u',0)]
	    self.lsModel.q[('grad(u)',0)] += self.massCorrModel.q[('grad(u)',0)]
	    self.lsModel.ebqe[('grad(u)',0)] += self.massCorrModel.ebqe[('grad(u)',0)]
            #vof
            if self.edgeBasedStabilizationMethods==False:
                self.massCorrModel.setMassQuadrature()
            #else setMassQuadratureEdgeBasedStabilizationMethods is called within specialized nolinear solver
            
            #self.vofModel.q[('u',0)] += self.massCorrModel.q[('r',0)]
            #####print "********************max VOF************************",max(self.vofModel.q[('u',0)].flat[:])
        if self.checkMass:
            logEvent("Phase 0 mass after mass correction (VOF) %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                                self.vofModel.q[('m',0)],
                                                                                                self.massCorrModel.mesh.nElements_owned),),level=2)
            logEvent("Phase 0 mass (primitive) after mass correction (LS) %12.5e" % (Norms.scalarSmoothedHeavisideDomainIntegral(self.epsFactHeaviside,
                                                                                                                self.massCorrModel.elementDiameter,
                                                                                                                self.vofModel.q['dV'],
                                                                                                                self.lsModel.q[('m',0)],
                                                                                                                self.massCorrModel.mesh.nElements_owned),),level=2)
            logEvent("Phase 0 mass (consistent) after mass correction (LS) %12.5e" % (self.massCorrModel.calculateMass(self.lsModel.q[('m',0)]),),level=2)
        copyInstructions = {}
        #get the waterline on the obstacle if option set in NCLS (boundary==7)
	self.lsModel.computeWaterline(t)
        return copyInstructions
    def evaluate(self,t,c):
        import math
        if c[('u',0)].shape == self.q_u_ls.shape:
            u_ls = self.q_u_ls
            H_vof = self.q_H_vof
        elif c[('u',0)].shape == self.ebqe_u_ls.shape:
            u_ls = self.ebqe_u_ls
            H_vof = self.ebqe_H_vof
        elif self.ebq_u_ls is not None and c[('u',0)].shape == self.ebq_u_ls.shape:
            u_ls = self.ebq_u_ls
            H_vof = self.ebq_H_vof
        else:
            #\todo trap errors in TransportCoefficients.py
            u_ls = None
            H_vof = None
        if u_ls is not None and H_vof is not None:
            if self.useC:
                if self.sd:
                    self.levelSetConservationCoefficientsEvaluate_sd(self.epsHeaviside,
                                                                     self.epsDirac,
                                                                     u_ls,
                                                                     H_vof,
                                                                     c[('u',0)],
                                                                     c[('r',0)],
                                                                     c[('dr',0,0)])
                else:
                    self.levelSetConservationCoefficientsEvaluate(self.epsHeaviside,
                                                                  self.epsDirac,
                                                                  self.epsDiffusion,
                                                                  u_ls,
                                                                  H_vof,
                                                                  c[('u',0)],
                                                                  c[('r',0)],
                                                                  c[('dr',0,0)],
                                                                  c[('a',0,0)])
        if (self.checkMass and c[('u',0)].shape == self.q_u_ls.shape):
            self.m_tmp[:] = H_vof
            self.m_tmp += self.massCorrModel.q[('r',0)]
            logEvent("mass correction during Newton %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                     self.massCorrModel.q[('r',0)],
                                                                                     self.massCorrModel.mesh.nElements_owned),),level=2)
            logEvent("Phase 0 mass during Newton %12.5e" % (Norms.scalarDomainIntegral(self.vofModel.q['dV'],
                                                                                 self.m_tmp,
                                                                                  self.massCorrModel.mesh.nElements_owned),),level=2)

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
                 bdyNullSpace=False):#,
        self.useConstantH=coefficients.useConstantH
        from proteus import Comm
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
        self.bdyNullSpace = bdyNullSpace
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
	if self.stabilization is not None:
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
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux is not None) or
                                                 (numericalFluxType is not None) or
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
        if self.stabilization is not None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
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
            #print self.nQuadraturePoints_element
            if self.nSpace_global == 3:
                assert(self.nQuadraturePoints_element == 5)
            elif self.nSpace_global == 2:
                assert(self.nQuadraturePoints_element == 6)
            elif self.nSpace_global == 1:
                assert(self.nQuadraturePoints_element == 3)

            #print self.nElementBoundaryQuadraturePoints_elementBoundary
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
        #uncomment this to store q arrays, see calculateElementQuadrature below 
        #self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('grad(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('r',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')

        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('grad(u)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')


        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
        logEvent(memory("element and element boundary Jacobians","OneLevelTransport"),level=4)
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

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        logEvent(memory("TimeIntegration","OneLevelTransport"),level=4)
        logEvent("Calculating numerical quadrature formulas",2)
        self.calculateQuadrature()
        self.setupFieldStrides()

        # (MQL)
        self.MassMatrix=None #consistent mass matrix
        self.LumpedMassMatrix=None
        self.rhs_mass_correction = None
        self.MassMatrix_sparseFactor=None
        self.Jacobian_sparseFactor=None
        self.lumped_L2p_vof_mass_correction=None
        self.limited_L2p_vof_mass_correction=None
        self.L2p_vof_mass_correction=None

        comm = Comm.get()
        self.comm=comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions,"You must use a numerical flux to apply weak boundary conditions for parallel runs"

        logEvent(memory("stride+offset","OneLevelTransport"),level=4)
        if numericalFluxType is not None:
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
        self.globalResidualDummy = None
        compKernelFlag = 0
        if self.coefficients.useConstantH:
            self.elementDiameter = self.mesh.elementDiametersArray.copy()
            self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
        else:
            self.elementDiameter = self.mesh.elementDiametersArray
        self.mcorr = cMCorr_base(self.nSpace_global,
                                 self.nQuadraturePoints_element,
                                 self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                 self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                 self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                                 compKernelFlag)
    #mwf these are getting called by redistancing classes,
    def FCTStep(self):
        rowptr, colind, MassMatrix = self.MassMatrix.getCSRrepresentation()
        if (self.limited_L2p_vof_mass_correction is None):
            self.limited_L2p_vof_mass_correction = numpy.zeros(self.LumpedMassMatrix.size,'d')

        self.mcorr.FCTStep(
            self.nnz, #number of non zero entries 
            len(rowptr)-1, #number of DOFs
            self.LumpedMassMatrix, #Lumped mass matrix
            self.L2p_vof_mass_correction, # high order projection
            self.lumped_L2p_vof_mass_correction, #low order projection
            self.limited_L2p_vof_mass_correction,
            rowptr, #Row indices for Sparsity Pattern (convenient for DOF loops)
            colind, #Column indices for Sparsity Pattern (convenient for DOF loops)
            MassMatrix)        
    
    def calculateCoefficients(self):
        pass
    def calculateElementResidual(self):
        if self.globalResidualDummy is not None:
            self.getResidual(self.u[0].dof,self.globalResidualDummy)
    def getResidual(self,u,r):
        import pdb
        import copy
        from proteus.flcbdfWrappers import globalSum
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.setUnknowns(u)


        #no flux boundary conditions
        self.mcorr.calculateResidual(#element
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
	    self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.coefficients.q_u_ls,
            self.coefficients.q_n_ls,
            self.coefficients.ebqe_u_ls,
            self.coefficients.ebqe_n_ls,
            self.coefficients.q_H_vof,
            self.q[('u',0)],
            self.q[('grad(u)',0)],
            self.ebqe[('u',0)],
            self.ebqe[('grad(u)',0)],
            self.q[('r',0)],
            self.coefficients.q_porosity,
            self.offset[0],self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray)
        logEvent("Global residual",level=9,data=r)
        self.coefficients.massConservationError = fabs(globalSum(r[:self.mesh.nNodes_owned].sum()))
        logEvent("   Mass Conservation Error: ",level=3,data=self.coefficients.massConservationError)
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy is None:
            self.globalResidualDummy = numpy.zeros(r.shape,'d')

    # GET MASS MATRIX # (MQL)
    def getMassMatrix(self):
        #NOTE. Both, the consistent and the lumped mass matrix must be init to zero
        if (self.MassMatrix is None):
            rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
            self.MassMatrix_a = nzval.copy()
            nnz = nzval.shape[-1] #number of non-zero entries in sparse matrix
            self.MassMatrix = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                           self.nFreeDOF_global[0],
                                                           nnz,
                                                           self.MassMatrix_a,
                                                           colind,
                                                           rowptr)
            # Lumped mass matrix
            self.LumpedMassMatrix = numpy.zeros(rowptr.size-1,'d')
        else:
            self.LumpedMassMatrix.fill(0.0)

        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       self.MassMatrix)
        self.mcorr.calculateMassMatrix(#element
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
            self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.coefficients.q_u_ls,
			self.coefficients.q_n_ls,
            self.coefficients.q_H_vof,
            self.coefficients.q_porosity,
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            self.MassMatrix, 
            self.LumpedMassMatrix)

    def getJacobian(self,jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,jacobian)
        self.mcorr.calculateJacobian(#element
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
            self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.coefficients.q_u_ls,
			self.coefficients.q_n_ls,
            self.coefficients.q_H_vof,
            self.coefficients.q_porosity,
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            jacobian)
        logEvent("Jacobian ",level=10,data=jacobian)
        #mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian

    def elementSolve(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.setUnknowns(u)

        #no flux boundary conditions
        self.mcorr.elementSolve(#element
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
	    self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.coefficients.q_u_ls,
            self.coefficients.q_n_ls,
            self.coefficients.ebqe_u_ls,
            self.coefficients.ebqe_n_ls,
            self.coefficients.q_H_vof,
            self.q[('u',0)],
            self.q[('grad(u)',0)],
            self.ebqe[('u',0)],
            self.ebqe[('grad(u)',0)],
            self.q[('r',0)],
            self.coefficients.q_porosity,
            self.offset[0],self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.maxIts,
            self.atol)

    def elementConstantSolve(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.setUnknowns(u)

        #no flux boundary conditions
        self.mcorr.elementConstantSolve(#element
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
	    self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.coefficients.q_u_ls,
            self.coefficients.q_n_ls,
            self.coefficients.ebqe_u_ls,
            self.coefficients.ebqe_n_ls,
            self.coefficients.q_H_vof,
            self.q[('u',0)],
            self.q[('grad(u)',0)],
            self.ebqe[('u',0)],
            self.ebqe[('grad(u)',0)],
            self.q[('r',0)],
            self.coefficients.q_porosity,
            self.offset[0],self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.maxIts,
            self.atol)

    def globalConstantRJ(self,u,r,U):
        from proteus.flcbdfWrappers import globalSum
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.setUnknowns(u)

        #no flux boundary conditions
        (R,J) = self.mcorr.globalConstantRJ(#element
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
            self.mesh.nElements_owned,
	    self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,
            self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            self.coefficients.q_u_ls,
            self.coefficients.q_n_ls,
            self.coefficients.ebqe_u_ls,
            self.coefficients.ebqe_n_ls,
            self.coefficients.q_H_vof,
            self.q[('u',0)],
            self.q[('grad(u)',0)],
            self.ebqe[('u',0)],
            self.ebqe[('grad(u)',0)],
            self.q[('r',0)],
            self.offset[0],self.stride[0],
            r,
            self.coefficients.q_porosity,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.maxIts,
            self.atol,
            U)
        R = globalSum(R)
        J = globalSum(J)
        self.coefficients.massConservationError = fabs(R)
        return (R,J)

    def globalConstantSolve(self,u,r):
        U=0.0
        R=0.0
        J=0.0
        (R,J) = self.globalConstantRJ(u,r,U)
        its=0
        logEvent("   Mass Conservation Residual 0 ",level=3,data=R)
        RNORM_OLD=fabs(R)
        while ((fabs(R) > self.atol and its < self.maxIts) or its<1):
            U -= R/(J+1.0e-8)
            (R,J) = self.globalConstantRJ(u,r,U)
            lsits=0
            while(fabs(R) > 0.99*RNORM_OLD and lsits < self.maxLSits):
                lsits +=1
                U += (0.5)**lsits * (R/(J+1.0e-8))
                (R,J) = self.globalConstantRJ(u,r,U)
            its+=1
            logEvent("   Mass Conservation Residual "+`its`+" ",level=3,data=R)
        self.u[0].dof.flat[:] = U
    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        #uncomment this to store q arrays
        #self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
        #                                          self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization is not None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        pass
    def calculateExteriorElementBoundaryQuadrature(self):
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
    def estimate_mt(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass
    def calculateMass(self,q_phi):
        from proteus.flcbdfWrappers import globalSum
        return globalSum(self.mcorr.calculateMass(#element
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
            self.mesh.nElements_owned,
	    self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.u[0].dof,
            q_phi,#self.coefficients.q_u_ls,
            self.coefficients.q_n_ls,
            self.coefficients.ebqe_u_ls,
            self.coefficients.ebqe_n_ls,
            self.coefficients.q_H_vof,
            self.q[('u',0)],
            self.q[('grad(u)',0)],
            self.ebqe[('u',0)],
            self.ebqe[('grad(u)',0)],
            self.q[('r',0)],
            self.coefficients.q_porosity,
            self.offset[0],self.stride[0],
            self.u[0].dof,#dummy r,not used
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray))

    def setMassQuadratureEdgeBasedStabilizationMethods(self):
        #Compute mass matrix
        #Set rhs of mass correction to zero
        if self.rhs_mass_correction is None:
            self.rhs_mass_correction = numpy.zeros(self.coefficients.vofModel.u[0].dof.shape,'d')
            self.lumped_L2p_vof_mass_correction = numpy.zeros(self.coefficients.vofModel.u[0].dof.shape,'d')
            self.L2p_vof_mass_correction = numpy.zeros(self.coefficients.vofModel.u[0].dof.shape,'d')
        else: 
            self.rhs_mass_correction.fill(0.0)
        
        self.mcorr.setMassQuadratureEdgeBasedStabilizationMethods(#element
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
	    self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.coefficients.lsModel.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.coefficients.lsModel.u[0].dof,
            self.coefficients.q_u_ls,
            self.coefficients.q_n_ls,
            self.coefficients.ebqe_u_ls,
            self.coefficients.ebqe_n_ls,
            self.coefficients.q_H_vof,
            self.q[('u',0)],
            self.q[('grad(u)',0)],
            self.ebqe[('u',0)],
            self.ebqe[('grad(u)',0)],
            self.q[('r',0)],
            self.coefficients.q_porosity,
            self.offset[0],self.stride[0],
            self.u[0].dof,#dummy r,not used
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.rhs_mass_correction, # (MQL): compute rhs for L2 projection. 
            self.lumped_L2p_vof_mass_correction, 
            self.LumpedMassMatrix,
            self.lumped_L2p_vof_mass_correction.size)
        
    def setMassQuadrature(self):
        self.mcorr.setMassQuadrature(#element
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
	    self.coefficients.useMetrics,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.coefficients.lsModel.u[0].femSpace.dofMap.l2g,
            self.elementDiameter,#self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.coefficients.lsModel.u[0].dof,
            self.coefficients.q_u_ls,
            self.coefficients.q_n_ls,
            self.coefficients.ebqe_u_ls,
            self.coefficients.ebqe_n_ls,
            self.coefficients.q_H_vof,
            self.q[('u',0)],
            self.q[('grad(u)',0)],
            self.ebqe[('u',0)],
            self.ebqe[('grad(u)',0)],
            self.q[('r',0)],
            self.coefficients.q_porosity,
            self.offset[0],self.stride[0],
            self.u[0].dof,#dummy r,not used
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.vofModel.u[0].dof)

    def calculateSolutionAtQuadrature(self):
        pass
    def updateAfterMeshMotion(self):
        pass


class DummyNewton(proteus.NonlinearSolvers.NonlinearSolver):
    def __init__(self,
                 linearSolver,
                 F,J=None,du=None,par_du=None,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True,
                 fullNewton=True,
                 directSolver=False,
                 EWtol=True,
                 maxLSits = 100):
        import copy
        self.par_du = par_du
        if par_du is not None:
            F.dim_proc = par_du.dim_proc
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.linearSolver = linearSolver
        self.directSolver = directSolver
        self.lineSearch = True
        #mwf turned back on self.lineSearch = False
        self.EWtol=EWtol
        #mwf added
        self.maxLSits = maxLSits
        if self.linearSolver.computeEigenvalues:
            self.JLast = copy.deepcopy(self.J)
            self.J_t_J = copy.deepcopy(self.J)
            self.dJ_t_dJ = copy.deepcopy(self.J)
            self.JLsolver=LU(self.J_t_J,computeEigenvalues=True)
            self.dJLsolver=LU(self.dJ_t_dJ,computeEigenvalues=True)
            self.u0 = numpy.zeros(self.F.dim,'d')
    def info(self):
        return "Not Implemented"
    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        self.F.q[('r',0)].flat[:]=0.0
        self.F.q[('u',0)].flat[:]=0.0
        self.failedFlag=False
        return self.failedFlag
class ElementNewton(proteus.NonlinearSolvers.NonlinearSolver):
    def __init__(self,
                 linearSolver,
                 F,J=None,du=None,par_du=None,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True,
                 fullNewton=True,
                 directSolver=False,
                 EWtol=True,
                 maxLSits = 100):
        import copy
        self.par_du = par_du
        if par_du is not None:
            F.dim_proc = par_du.dim_proc
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.linearSolver = linearSolver
        self.directSolver = directSolver
        self.lineSearch = True
        #mwf turned back on self.lineSearch = False
        self.EWtol=EWtol
        #mwf added
        self.maxLSits = maxLSits
        if self.linearSolver.computeEigenvalues:
            self.JLast = copy.deepcopy(self.J)
            self.J_t_J = copy.deepcopy(self.J)
            self.dJ_t_dJ = copy.deepcopy(self.J)
            self.JLsolver=LU(self.J_t_J,computeEigenvalues=True)
            self.dJLsolver=LU(self.dJ_t_dJ,computeEigenvalues=True)
            self.u0 = numpy.zeros(self.F.dim,'d')
    def info(self):
        return "Not Implemented"
    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        self.F.maxIts = self.maxIts
        self.F.maxLSits = self.maxLSits
        self.F.atol = self.atol_r
        self.F.elementSolve(u,r)
        self.failedFlag=False
        return self.failedFlag
class ElementConstantNewton(proteus.NonlinearSolvers.NonlinearSolver):
    def __init__(self,
                 linearSolver,
                 F,J=None,du=None,par_du=None,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True,
                 fullNewton=True,
                 directSolver=False,
                 EWtol=True,
                 maxLSits = 100):
        import copy
        self.par_du = par_du
        if par_du is not None:
            F.dim_proc = par_du.dim_proc
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.linearSolver = linearSolver
        self.directSolver = directSolver
        self.lineSearch = True
        #mwf turned back on self.lineSearch = False
        self.EWtol=EWtol
        #mwf added
        self.maxLSits = maxLSits
        if self.linearSolver.computeEigenvalues:
            self.JLast = copy.deepcopy(self.J)
            self.J_t_J = copy.deepcopy(self.J)
            self.dJ_t_dJ = copy.deepcopy(self.J)
            self.JLsolver=LU(self.J_t_J,computeEigenvalues=True)
            self.dJLsolver=LU(self.dJ_t_dJ,computeEigenvalues=True)
            self.u0 = numpy.zeros(self.F.dim,'d')
    def info(self):
        return "Not Implemented"
    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        self.F.maxIts = self.maxIts
        self.F.maxLSits = self.maxLSits
        self.F.atol = self.atol_r
        self.F.elementConstantSolve(u,r)
        self.failedFlag=False
        return self.failedFlag
class GlobalConstantNewton(proteus.NonlinearSolvers.NonlinearSolver):
    def __init__(self,
                 linearSolver,
                 F,J=None,du=None,par_du=None,
                 rtol_r  = 1.0e-4,
                 atol_r  = 1.0e-16,
                 rtol_du = 1.0e-4,
                 atol_du = 1.0e-16,
                 maxIts  = 100,
                 norm = l2Norm,
                 convergenceTest = 'r',
                 computeRates = True,
                 printInfo = True,
                 fullNewton=True,
                 directSolver=False,
                 EWtol=True,
                 maxLSits = 100):
        import copy
        self.par_du = par_du
        if par_du is not None:
            F.dim_proc = par_du.dim_proc
        NonlinearSolver.__init__(self,F,J,du,
                                 rtol_r,
                                 atol_r,
                                 rtol_du,
                                 atol_du,
                                 maxIts,
                                 norm,
                                 convergenceTest,
                                 computeRates,
                                 printInfo)
        self.updateJacobian=True
        self.fullNewton=fullNewton
        self.linearSolver = linearSolver
        self.directSolver = directSolver
        self.lineSearch = True
        #mwf turned back on self.lineSearch = False
        self.EWtol=EWtol
        #mwf added
        self.maxLSits = maxLSits
        if self.linearSolver.computeEigenvalues:
            self.JLast = copy.deepcopy(self.J)
            self.J_t_J = copy.deepcopy(self.J)
            self.dJ_t_dJ = copy.deepcopy(self.J)
            self.JLsolver=LU(self.J_t_J,computeEigenvalues=True)
            self.dJLsolver=LU(self.dJ_t_dJ,computeEigenvalues=True)
            self.u0 = numpy.zeros(self.F.dim,'d')
    def info(self):
        return "Not Implemented"
    def solve(self,u,r=None,b=None,par_u=None,par_r=None):
        self.F.maxIts = self.maxIts
        self.F.maxLSits = self.maxLSits
        self.F.atol = self.atol_r
        self.F.globalConstantSolve(u,r)
        self.failedFlag=False
        return self.failedFlag

from proteus.flcbdfWrappers import globalSum

def conservationNorm(x):
    return fabs(globalSum(sum(x.flat)))

class Newton_controller(proteus.StepControl.Newton_controller):
    def __init__(self,model,nOptions):
        proteus.StepControl.Newton_controller.__init__(self,model,nOptions)
    def initializeTimeHistory(self):
        proteus.StepControl.Newton_controller.initializeTimeHistory(self)
        for m,u,r in zip(self.model.levelModelList,
                         self.model.uList,
                         self.model.rList):
            u.flat[:]=0.0
            m.getResidual(u,r)
            m.coefficients.postStep(self.t_model)
            m.coefficients.vofModel.updateTimeHistory(self.t_model,resetFromDOF=False)
            m.coefficients.vofModel.timeIntegration.updateTimeHistory(resetFromDOF=False)
