from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import proteus
from proteus.mprans.cCLSVOF import *

class NumericalFlux(proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                                        getAdvectiveFluxBoundaryConditions,
                                                                                        getDiffusiveFluxBoundaryConditions)

class Coefficients(proteus.TransportCoefficients.TC_base):
    def __init__(self,
                 LS_model=None,
                 V_model=0,
                 RD_model=None,
                 ME_model=1,
                 VOS_model=None,
                 checkMass=False,
                 epsFact=0.0,
                 useMetrics=0.0,
                 setParamsFunc=None,
                 movingDomain=False,
                 forceStrongConditions=0,
                 # OUTPUT quantDOFs
                 outputQuantDOFs = True, # mql. I use it to visualize H(u) at the DOFs
                 computeMetrics = 0, #0, 1, 2 or 3
                 # SPIN UP STEP #
                 doSpinUpStep=False, # To achieve high order with Bernstein polynomials
                 disc_ICs=False, # Is the init condition a characteristic function?
                 # NONLINEAR CLSVOF
                 timeOrder=1,
                 epsFactHeaviside=1.5,
                 epsFactDirac=1.5,
                 epsFactRedist=0.33,
                 lambdaFact=1.0,
                 alpha='inf'): #lambda parameter in CLSVOF paper
        assert timeOrder==1, "timeOrder must be 1. It will be deleted after 1st paper"
        assert timeOrder==1 or timeOrder==2, "timeOrder must be 1 or 2"
        assert computeMetrics in [0,1,2,3]
        # 0: don't compute metrics
        # 1: compute change in volume at ETS (every time step)
        # 2: compute several metrics at ETS (every time step)
        # 3: compute metrics at EOS (end of simulations). Needs an exact solution
        self.useMetrics=useMetrics
        self.doSpinUpStep=doSpinUpStep
        self.disc_ICs=disc_ICs
        self.timeOrder=timeOrder
        self.computeMetrics=computeMetrics
        self.epsFactHeaviside=epsFactHeaviside
        self.epsFactDirac=epsFactDirac
        self.epsFactRedist=epsFactRedist
        self.lambdaFact=lambdaFact
        self.variableNames=['clsvof']
        nc=1
        mass={0:{0:'linear'}}
        advection={0:{0:'linear'}}
        hamiltonian={}
        diffusion={}
        potential={}
        reaction={}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         self.variableNames,
                         movingDomain=movingDomain)
        self.epsFact=epsFact
        self.flowModelIndex=V_model
        self.modelIndex=ME_model
        self.RD_modelIndex=RD_model
        self.VOS_model=VOS_model
        self.checkMass = checkMass
        #VRANS
        self.setParamsFunc   = setParamsFunc
        self.flowCoefficients=None
        self.movingDomain=movingDomain
        self.forceStrongConditions=forceStrongConditions
        self.outputQuantDOFs=outputQuantDOFs
        self.alpha=alpha
        self.freeze_interface_during_preRedistancing = False
        if self.alpha=='inf':
            self.alpha = 0
            self.freeze_interface_during_preRedistancing = True
        # VRANS
        self.q_vos = None
        self.ebqe_vos = None

    def initializeMesh(self,mesh):
        self.eps = self.epsFact*mesh.h
    def attachModels(self,modelList):
        #self
        self.model = modelList[self.modelIndex]
        #flow model
        if self.flowModelIndex is not None:
            self.flowCoefficients = modelList[self.flowModelIndex].coefficients
            if ('velocity',0) in modelList[self.flowModelIndex].q:
                self.q_v = modelList[self.flowModelIndex].q[('velocity',0)]
                self.ebqe_v = modelList[self.flowModelIndex].ebqe[('velocity',0)]
            else:
                self.q_v = modelList[self.flowModelIndex].q[('f',0)]
                self.ebqe_v = modelList[self.flowModelIndex].ebqe[('f',0)]
            if ('velocity',0) in modelList[self.flowModelIndex].ebq:
                self.ebq_v = modelList[self.flowModelIndex].ebq[('velocity',0)]
            else:
                if ('f',0) in modelList[self.flowModelIndex].ebq:
                    self.ebq_v = modelList[self.flowModelIndex].ebq[('f',0)]
        self.q_v_old = numpy.copy(self.q_v)
        self.q_v_tStar = numpy.copy(self.q_v)
        #VRANS
        if self.VOS_model is not None:
            self.model.q_vos = modelList[self.VOS_model].q[('u',0)]
            self.model.ebqe_vos = modelList[self.VOS_model].ebqe[('u',0)]
            self.q_vos = self.model.q_vos
            self.ebqe_vos = self.model.ebqe_vos
        else:
            q_porosity = numpy.ones(modelList[self.modelIndex].q[('u',0)].shape,'d')
            ebqe_porosity = numpy.ones(self.model.ebqe[('u', 0)].shape, 'd')
            # set porosity from setParamsFunc
            if self.setParamsFunc is not None:
                self.setParamsFunc(modelList[self.modelIndex].q['x'],q_porosity)
                self.setParamsFunc(modelList[self.modelIndex].ebqe['x'],ebqe_porosity)
            # set porosity from flowCoefficients
            elif self.flowCoefficients is not None:
                if hasattr(self.flowCoefficients,'q_porosity'):
                    q_porosity[:] = self.flowCoefficients.q_porosity
                if hasattr(self.flowCoefficients,'ebqe_porosity'):
                    ebqe_porosity[:] = self.flowCoefficients.ebqe_porosity
            #
            self.q_vos[:] = 1. - q_porosity[:]
            self.ebqe_vos[:] = 1. - ebqe_porosity[:]

    def initializeElementQuadrature(self,t,cq):
        if self.flowModelIndex == None:
            self.q_v = numpy.ones(cq[('f',0)].shape,'d')
        #VRANS
        self.q_vos = numpy.zeros(cq[('u',0)].shape,'d')

    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        if self.flowModelIndex == None:
            self.ebq_v = numpy.ones(cebq[('f',0)].shape,'d')
        #VRANS
        self.ebq_vos = numpy.zeros(cebq[('u',0)].shape,'d')

    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        if self.flowModelIndex == None:
            self.ebqe_v = numpy.ones(cebqe[('f',0)].shape,'d')
        #VRANS
        self.ebqe_vos = numpy.zeros(cebqe[('u',0)].shape,'d')

    def preStep(self,t,firstStep=False):
        self.model.getResidualBeforeFirstStep = False
        # SAVE INITIAL CONDITION TO MEASURE ERRORS #
        if (self.computeMetrics > 0 and firstStep==True):
            # Overwritten if spin up step is taken
            self.model.u0_dof[:] = self.model.u[0].dof
        # COMPUTE NEW VELOCITY (if given by user) #
        if self.model.hasVelocityFieldAsFunction:
            self.model.updateVelocityFieldAsFunction()
        if (firstStep==True):
            self.q_v_old[:] = self.q_v
        # GET MAX VELOCITY #
        self.VelMax = max(self.q_v.max(),1E-6)
        copyInstructions = {}
        return copyInstructions

    def postStep(self,t,firstStep=False):
        self.model.quantDOFs[:] = self.model.interface_locator
        # SAVE OLD SOLUTION #
        self.model.u_dof_old[:] = self.model.u[0].dof
        # SAVE OLD VELOCITY #
        self.q_v_old[:] = self.q_v
        # NORM FACTOR #
        self.model.norm_factor_lagged=np.maximum(self.model.max_distance - self.model.mean_distance,
                                                 self.model.mean_distance - self.model.min_distance)

        # Compute metrics at end of time step
        if self.computeMetrics == 1 or self.computeMetrics == 2: #compute metrics at ETS
            self.model.getMetricsAtETS()
            if self.model.comm.isMaster():
                if self.computeMetrics == 1:
                    self.model.metricsAtETS.write(repr(self.model.timeIntegration.t)[:4]+",\t"+
                                                  repr(self.model.global_sV_err)+
                                                  "\n")
                    self.model.metricsAtETS.flush()
                else:
                    self.model.metricsAtETS.write(repr(self.model.timeIntegration.t)[:4]+","+
                                                  repr(self.model.timeIntegration.dt)+","+
                                                  repr(self.model.newton_iterations_stage1)+","+
                                                  repr(self.model.newton_iterations_stage2)+","+
                                                  repr(math.sqrt(self.model.global_R))+","+
                                                  repr(math.sqrt(self.model.global_sR))+","+
                                                  repr(self.model.global_V_err)+","+
                                                  repr(self.model.global_sV_err)+","+
                                                  repr(self.model.global_D_err)+
                                                  "\n")
                    self.model.metricsAtETS.flush()
        #
        self.model.q['dV_last'][:] = self.model.q['dV']
        copyInstructions = {}
        return copyInstructions

    def updateToMovingDomain(self,t,c):
        #in a moving domain simulation the velocity coming in is already for the moving domain
        pass

    def evaluate(self,t,c):
        pass

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

        self.mass_history=[]
        #
        #set the objects describing the method and boundary conditions
        #
        self.bdyNullSpace=bdyNullSpace
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
                if ci in coefficients.mass:
                    for flag in list(coefficients.mass[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if  ci in coefficients.advection:
                    for  flag  in list(coefficients.advection[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if  ci in coefficients.diffusion:
                    for diffusionDict in list(coefficients.diffusion[ci].values()):
                        for  flag  in list(diffusionDict.values()):
                            if flag != 'constant':
                                self.stabilizationIsNonlinear=True
                if  ci in coefficients.potential:
                    for flag in list(coefficients.potential[ci].values()):
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if ci in coefficients.reaction:
                    for flag in list(coefficients.reaction[ci].values()):
                        if  flag == 'nonlinear':
                            self.stabilizationIsNonlinear=True
                if ci in coefficients.hamiltonian:
                    for flag in list(coefficients.hamiltonian[ci].values()):
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
        self.nDOF_trial_element     = [u_j.femSpace.max_nDOF_element for  u_j in list(self.u.values())]
        self.nDOF_phi_trial_element     = [phi_k.femSpace.max_nDOF_element for  phi_k in list(self.phi.values())]
        self.n_phi_ip_element = [phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for  phi_k in list(self.phi.values())]
        self.nDOF_test_element     = [femSpace.max_nDOF_element for femSpace in list(self.testSpace.values())]
        self.nFreeDOF_global  = [dc.nFreeDOF_global for dc in list(self.dirichletConditions.values())]
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
                if I in elementQuadrature:
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization != None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if I in elementQuadrature:
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',)+I[1:]] = elementQuadrature
        if self.shockCapturing != None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if ('numDiff',ci,ci) in elementQuadrature:
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature[('numDiff',ci,ci)]
                    else:
                        elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff',ci,ci)] = elementQuadrature
        if massLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('m',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        if reactionLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('r',ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global,1)
        elementBoundaryQuadratureDict={}
        if isinstance(elementBoundaryQuadrature,dict): #set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if I in elementBoundaryQuadrature:
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
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('H(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('mH(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dV_u',0)] = (old_div(1.0,self.mesh.nElements_global))*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('grad(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('m',0)] = self.q[('u',0)]
        self.q[('m_last',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('mt',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['dV'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['dV_last'] = -1000*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',0)] = self.q[('u',0)].copy()
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',0,0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('H(u)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('grad(u)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')

        # JACOBIANS (FOR ELEMENT TRANSFORMATION)
        self.q[('J')] = np.zeros((self.mesh.nElements_global,
                                  self.nQuadraturePoints_element,
                                  self.nSpace_global,
                                  self.nSpace_global),
                                 'd')
        self.q[('inverse(J)')] = np.zeros((self.mesh.nElements_global,
                                           self.nQuadraturePoints_element,
                                           self.nSpace_global,
                                           self.nSpace_global),
                                          'd')
        self.q[('det(J)')] = np.zeros((self.mesh.nElements_global,
                                       self.nQuadraturePoints_element),
                                      'd')
        self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                         self.q['J'],
                                                         self.q['inverse(J)'],
                                                         self.q['det(J)'])
        self.q['abs(det(J))'] = np.abs(self.q['det(J)'])

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
        self.setupFieldStrides()

        #cek adding empty data member for low order numerical viscosity structures here for now
        self.Jacobian_sparseFactor=None

        comm = Comm.get()
        self.comm=comm
        if comm.size() > 1:
            assert numericalFluxType != None and numericalFluxType.useWeakDirichletConditions,"You must use a numerical flux to apply weak boundary conditions for parallel runs"

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
        if 'penalty' in self.ebq_global:
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN,k] = old_div(self.numericalFlux.penalty_constant,(self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power))
        #penalty term
        #cek move  to Numerical flux initialization
        if 'penalty' in self.ebqe:
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE,k] = old_div(self.numericalFlux.penalty_constant,self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
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
        for ci,fbcObject  in list(self.fluxBoundaryConditionsObjectsDict.items()):
            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc',ci)].shape,'i')
            for t,g in list(fbcObject.advectiveFluxBoundaryConditionsDict.items()):
                if ci in self.coefficients.advection:
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
        self.clsvof = cCLSVOF_base(self.nSpace_global,
                             self.nQuadraturePoints_element,
                             self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                             self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                             self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                             compKernelFlag)
        # strong Dirichlet boundary conditions
        self.forceStrongConditions=False
        if self.forceStrongConditions:
            self.dirichletConditionsForceDOF = DOFBoundaryConditions(self.u[0].femSpace,dofBoundaryConditionsSetterDict[0],weakDirichletConditions=False)

        if self.movingDomain:
            self.MOVING_DOMAIN=1.0
        else:
            self.MOVING_DOMAIN=0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape,'d')

        ################
        # SOME ASSERTS #
        ################
        assert isinstance(self.timeIntegration,proteus.TimeIntegration.BackwardEuler_cfl), "Use BackwardEuler_cfl"
        assert options.levelNonlinearSolver == proteus.NonlinearSolvers.CLSVOFNewton, "Use levelNonlinearSolver=CLSVOFNewton"

        #################
        # GENERAL STUFF #
        #################
        self.u_dof_old = None
        # Display info about CFL
        self.displayCFL=False
        self.timeStage=1 # stage of time integration

        ########################
        # POROSITY AS FUNCTION #
        ########################
        if ('porosityFieldAsFunction') in dir (options):
            X = {0:self.q[('x')][:,:,0],
                 1:self.q[('x')][:,:,1],
                 2:self.q[('x')][:,:,2]}
            self.coefficients.q_vos[:] = 1.-options.porosityFieldAsFunction(X)

            # BOUNDARY
            ebqe_X = {0:self.ebqe['x'][:,:,0],
                      1:self.ebqe['x'][:,:,1],
                      2:self.ebqe['x'][:,:,2]}
            self.coefficients.ebqe_vos[:] = 1.-options.porosityFieldAsFunction(ebqe_X)

        ################################
        # VELOCITY FIELD AS A FUNCTION #
        ################################
        self.getResidualBeforeFirstStep = True
        self.hasVelocityFieldAsFunction = False
        if ('velocityFieldAsFunction') in dir (options):
            self.velocityFieldAsFunction = options.velocityFieldAsFunction
            self.hasVelocityFieldAsFunction = True

        ##############################
        # VISUALIZATION OF VOF FIELD #
        ##############################
        self.vofDOFs = numpy.zeros(self.u[0].dof.shape,'d')
        self.par_vofDOFs = None
        self.lumped_mass_matrix = None
        # Aux quantity at DOFs
        self.quantDOFs = numpy.zeros(self.u[0].dof.shape,'d')

        #############################
        # L2 PROJECTION OF SOLUTION #
        #############################
        self.rhs_l2_proj = numpy.zeros(self.u[0].dof.shape,'d')
        self.projected_disc_ICs = numpy.zeros(self.u[0].dof.shape,'d')
        self.par_projected_disc_ICs = None

        from proteus.Comm import globalMax
        self.he_for_disc_ICs = 0.5*(-globalMax(-self.mesh.elementDiametersArray.min()) +
                                    globalMax(self.mesh.elementDiametersArray.max()))
        ###################################
        # PROJECTED NORMAL RECONSTRUCTION #
        ###################################
        self.consistentNormalReconstruction=False # for now it will always be False
        self.degree_polynomial=1
        try:
            self.degree_polynomial = self.u[0].femSpace.order
        except:
            pass
        #if self.degree_polynomial>1:
        #    self.consistentNormalReconstruction=True
        # For (lambda) normalization factor
        self.min_distance = 0.
        self.max_distance = 0.
        self.mean_distance = 0.
        self.volume_domain = 1.
        self.norm_factor_lagged = 1.0
        self.VelMax = 1.0
        self.weighted_lumped_mass_matrix = numpy.zeros(self.u[0].dof.shape,'d')
        # rhs for normal reconstruction
        self.rhs_qx = numpy.zeros(self.u[0].dof.shape,'d')
        self.rhs_qy = numpy.zeros(self.u[0].dof.shape,'d')
        self.rhs_qz = numpy.zeros(self.u[0].dof.shape,'d')
        # normal reconstruction
        self.projected_qx_tn = numpy.zeros(self.u[0].dof.shape,'d')
        self.projected_qy_tn = numpy.zeros(self.u[0].dof.shape,'d')
        self.projected_qz_tn = numpy.zeros(self.u[0].dof.shape,'d')
        self.projected_qx_tStar = numpy.zeros(self.u[0].dof.shape,'d')
        self.projected_qy_tStar = numpy.zeros(self.u[0].dof.shape,'d')
        self.projected_qz_tStar = numpy.zeros(self.u[0].dof.shape,'d')
        # parallel vectors for normal reconstruction
        self.par_projected_qx_tn = None
        self.par_projected_qy_tn = None
        self.par_projected_qz_tn = None
        self.par_projected_qx_tStar = None
        self.par_projected_qy_tStar = None
        self.par_projected_qz_tStar = None
        # Interface locator
        self.interface_locator = numpy.zeros(self.u[0].dof.shape,'d')
        self.preRedistancingStage = 0

        ###########################
        # CREATE PARALLEL VECTORS #
        ###########################
        #n=self.mesh.subdomainMesh.nNodes_owned
        #N=self.mesh.nNodes_global
        #nghosts=self.mesh.subdomainMesh.nNodes_global - self.mesh.subdomainMesh.nNodes_owned
        n=self.u[0].par_dof.dim_proc
        N=self.u[0].femSpace.dofMap.nDOF_all_processes
        nghosts = self.u[0].par_dof.nghosts
        subdomain2global=self.u[0].femSpace.dofMap.subdomain2global
        self.par_projected_qx_tn = proteus.LinearAlgebraTools.ParVec_petsc4py(self.projected_qx_tn,
                                                                              bs=1,
                                                                              n=n,N=N,nghosts=nghosts,
                                                                              subdomain2global=subdomain2global)
        self.par_projected_qy_tn = proteus.LinearAlgebraTools.ParVec_petsc4py(self.projected_qy_tn,
                                                                              bs=1,
                                                                              n=n,N=N,nghosts=nghosts,
                                                                              subdomain2global=subdomain2global)
        self.par_projected_qz_tn = proteus.LinearAlgebraTools.ParVec_petsc4py(self.projected_qz_tn,
                                                                              bs=1,
                                                                              n=n,N=N,nghosts=nghosts,
                                                                              subdomain2global=subdomain2global)
        self.par_projected_qx_tStar = proteus.LinearAlgebraTools.ParVec_petsc4py(self.projected_qx_tStar,
                                                                                 bs=1,
                                                                                 n=n,N=N,nghosts=nghosts,
                                                                                 subdomain2global=subdomain2global)
        self.par_projected_qy_tStar = proteus.LinearAlgebraTools.ParVec_petsc4py(self.projected_qy_tStar,
                                                                                 bs=1,
                                                                                 n=n,N=N,nghosts=nghosts,
                                                                                 subdomain2global=subdomain2global)
        self.par_projected_qz_tStar = proteus.LinearAlgebraTools.ParVec_petsc4py(self.projected_qz_tStar,
                                                                                 bs=1,
                                                                                 n=n,N=N,nghosts=nghosts,
                                                                                 subdomain2global=subdomain2global)
        #
        self.par_vofDOFs = proteus.LinearAlgebraTools.ParVec_petsc4py(self.vofDOFs,
                                                                      bs=1,
                                                                      n=n,N=N,nghosts=nghosts,
                                                                      subdomain2global=subdomain2global)
        self.par_projected_disc_ICs = proteus.LinearAlgebraTools.ParVec_petsc4py(self.projected_disc_ICs,
                                                                                 bs=1,
                                                                                 n=n,N=N,nghosts=nghosts,
                                                                                 subdomain2global=subdomain2global)

        ################
        # SPIN UP STEP #
        ################
        self.spinUpStepTaken=False
        self.uInitial = None
        if self.coefficients.doSpinUpStep:
            self.uInitial = numpy.zeros(self.q[('u',0)].shape,'d')
            X = {0:self.q[('x')][:,:,0],
                 1:self.q[('x')][:,:,1],
                 2:self.q[('x')][:,:,2]}
            self.uInitial[:] = options.initialConditions[0].uOfXT(X,0)

        ###########
        # METRICS #
        ###########
        self.hasExactSolution = False
        if ('exactSolution') in dir (options):
            self.hasExactSolution = True
            self.exactSolution = options.exactSolution

        self.u0_dof = numpy.copy(self.u[0].dof)
        self.newton_iterations_stage1 = 0.0
        self.newton_iterations_stage2 = 0.0
        # for interface quality
        self.global_I_err = 0.0
        self.global_sI_err = 0.0
        # for residual of conservation law
        self.R_vector = numpy.zeros(self.u[0].dof.shape,'d')
        self.sR_vector = numpy.zeros(self.u[0].dof.shape,'d')
        self.global_R = 0.0
        self.global_sR = 0.0
        # for conservation of volume
        self.global_V = 0.0
        self.global_V0 = 0.0
        self.global_sV = 0.0
        self.global_sV0 = 0.0
        self.global_V_err = 0.0
        self.global_sV_err = 0.0
        # for distance property
        self.global_D_err = 0.0
        self.global_L2_err = 0.0
        self.global_L2Banded_err = 0.0
        self.global_sH_L2_err = 0.0
        if self.coefficients.computeMetrics > 0 and self.comm.isMaster():
            if self.hasExactSolution and self.coefficients.computeMetrics==3: # at EOS
                self.metricsAtEOS = open(self.name+"_metricsAtEOS.csv","w")
                self.metricsAtEOS.write('global_I_err'+","+
                                        'global_sI_err'+","+
                                        'global_V_err'+","+
                                        'global_sV_err'+","+
                                        'global_D_err'+","+
                                        'global_L2_err'+","+
                                        'global_L2Banded_err'+","+
                                        'global_sH_L2_err'+"\n")
            elif self.coefficients.computeMetrics in [1,2]:
                self.metricsAtETS = open(self.name+"_metricsAtETS.csv","w")
                if self.coefficients.computeMetrics==1:
                    self.metricsAtETS.write('time'+","+
                                            'global_sV_err'+
                                            "\n")
                else:
                    self.metricsAtETS.write('time'+","+
                                            'time_step'+","+
                                            'newton_iterations_stage1'+","+
                                            'newton_iterations_stage2'+","+
                                            'global_R'+","+
                                            'global_sR'+","+
                                            'global_V_err'+","+
                                            'global_sV_err'+","+
                                            'global_D_err'+
                                            "\n")

    #mwf these are getting called by redistancing classes,
    def calculateCoefficients(self):
        pass

    def assembleSpinUpSystem(self,residual,jacobian):
        residual.fill(0.0)
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        self.clsvof.assembleSpinUpSystem(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.psi,
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.uInitial,
            self.offset[0],self.stride[0],
            residual,
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            jacobian)

    def updateVelocityFieldAsFunction(self):
        X = {0:self.q[('x')][:,:,0],
             1:self.q[('x')][:,:,1],
             2:self.q[('x')][:,:,2]}
        time = self.timeIntegration.t
        self.coefficients.q_v[...,0] = self.velocityFieldAsFunction[0](X,time)
        self.coefficients.q_v[...,1] = self.velocityFieldAsFunction[1](X,time)
        if (self.nSpace_global==3):
            self.coefficients.q_v[...,2] = self.velocityFieldAsFunction[2](X,time)

        # BOUNDARY
        ebqe_X = {0:self.ebqe['x'][:,:,0],
                  1:self.ebqe['x'][:,:,1],
                  2:self.ebqe['x'][:,:,2]}
        self.coefficients.ebqe_v[...,0] = self.velocityFieldAsFunction[0](ebqe_X,time)
        self.coefficients.ebqe_v[...,1] = self.velocityFieldAsFunction[1](ebqe_X,time)
        if (self.nSpace_global==3):
            self.coefficients.ebqe_v[...,2] = self.velocityFieldAsFunction[2](ebqe_X,time)

    def runAtEOS(self):
        if self.coefficients.computeMetrics ==3 and self.hasExactSolution:
            # Get exact solution at quad points
            u_exact = numpy.zeros(self.q[('u',0)].shape,'d')
            X = {0:self.q[('x')][:,:,0],
                 1:self.q[('x')][:,:,1],
                 2:self.q[('x')][:,:,2]}
            t = self.timeIntegration.t
            u_exact[:] = self.exactSolution[0](X,t)
            self.getMetricsAtEOS(u_exact)
            if self.comm.isMaster():
                self.metricsAtEOS.write(repr(self.global_I_err)+","+
                                        repr(self.global_sI_err)+","+
                                        repr(self.global_V_err)+","+
                                        repr(self.global_sV_err)+","+
                                        repr(self.global_D_err)+","+
                                        repr(np.sqrt(self.global_L2_err))+","+
                                        repr(np.sqrt(self.global_L2Banded_err))+","+
                                        repr(np.sqrt(self.global_sH_L2_err))+"\n")
                self.metricsAtEOS.flush()

    ####################################3
    def getMetricsAtETS(self): #ETS=Every Time Step
        """
        Calculate some metrics at ETS (Every Time Step)
        """

        self.R_vector.fill(0.)
        self.sR_vector.fill(0.)
        (global_V,
         global_V0,
         global_sV,
         global_sV0,
         global_D_err) = self.clsvof.calculateMetricsAtETS(
             self.timeIntegration.dt,
             self.u[0].femSpace.elementMaps.psi,
             self.u[0].femSpace.elementMaps.grad_psi,
             self.mesh.nodeArray,
             self.mesh.elementNodesArray,
             self.elementQuadratureWeights[('u',0)],
             self.u[0].femSpace.psi,
             self.u[0].femSpace.grad_psi,
             self.u[0].femSpace.psi,
             #physics
             self.mesh.nElements_global,
             self.mesh.nElements_owned,
             self.coefficients.useMetrics,
             self.coefficients.q_vos,
             self.u[0].femSpace.dofMap.l2g,
             self.mesh.elementDiametersArray,
             self.mesh.nodeDiametersArray,
             self.degree_polynomial,
             self.coefficients.epsFactHeaviside,
             self.u[0].dof, # This unp1
             self.u_dof_old,
             self.u0_dof,
             self.coefficients.q_v,
             self.offset[0],self.stride[0],
             self.nFreeDOF_global[0], #numDOFs
             self.R_vector,
             self.sR_vector)

        from proteus.Comm import globalSum
        # metrics about conservation
        self.global_V = globalSum(global_V)
        self.global_V0 = globalSum(global_V0)
        self.global_sV = globalSum(global_sV)
        self.global_sV0 = globalSum(global_sV0)
        self.global_V_err = old_div(np.abs(self.global_V-self.global_V0),self.global_V0)
        self.global_sV_err = old_div(np.abs(self.global_sV-self.global_sV0),self.global_sV0)
        # metrics about distance property
        self.global_D_err = globalSum(global_D_err)
        # compute global_R and global_sR
        n=self.mesh.subdomainMesh.nNodes_owned
        self.global_R = np.sqrt(globalSum(np.dot(self.R_vector[0:n],self.R_vector[0:n])))
        self.global_sR = np.sqrt(globalSum(np.dot(self.sR_vector[0:n],self.sR_vector[0:n])))

    def getMetricsAtEOS(self,u_exact): #EOS=End Of Simulation
        import copy
        """
        Calculate some metrics at EOS (End Of Simulation)
        """

        (global_I_err,
         global_sI_err,
         global_V,
         global_V0,
         global_sV,
         global_sV0,
         global_D_err,
         global_L2_err,
         global_L2Banded_err,
         global_area_band,
         global_sH_L2_err) = self.clsvof.calculateMetricsAtEOS(#element
             self.u[0].femSpace.elementMaps.psi,
             self.u[0].femSpace.elementMaps.grad_psi,
             self.mesh.nodeArray,
             self.mesh.elementNodesArray,
             self.elementQuadratureWeights[('u',0)],
             self.u[0].femSpace.psi,
             self.u[0].femSpace.grad_psi,
             self.u[0].femSpace.psi,
             #physics
             self.mesh.nElements_global,
             self.mesh.nElements_owned,
             self.coefficients.useMetrics,
             self.u[0].femSpace.dofMap.l2g,
             self.mesh.elementDiametersArray,
             self.mesh.nodeDiametersArray,
             self.degree_polynomial,
             self.coefficients.epsFactHeaviside,
             self.u[0].dof, # This is u_lstage due to update stages in RKEV
             self.u0_dof,
             u_exact,
             self.offset[0],self.stride[0])

        from proteus.Comm import globalSum
        # Interface metrics
        self.global_I_err = globalSum(global_I_err)
        self.global_sI_err = globalSum(global_sI_err)
        # conservation metrics
        self.global_V = globalSum(global_V)
        self.global_V0 = globalSum(global_V0)
        self.global_sV = globalSum(global_sV)
        self.global_sV0 = globalSum(global_sV0)
        self.global_V_err = old_div(np.abs(self.global_V-self.global_V0),self.global_V0)
        self.global_sV_err = old_div(np.abs(self.global_sV-self.global_sV0),self.global_sV0)
        # distance property metric
        self.global_D_err = globalSum(global_D_err)
        # L2 error on level set
        self.global_L2_err = globalSum(global_L2_err)
        self.global_L2Banded_err = old_div(globalSum(global_L2Banded_err),globalSum(global_area_band))
        self.global_sH_L2_err = globalSum(global_sH_L2_err)

    ###############################################

    def calculateElementResidual(self):
        if self.globalResidualDummy != None:
            self.getResidual(self.u[0].dof,self.globalResidualDummy)

    def FCTStep(self,
                limited_solution,
                soln,
                low_order_solution,
                high_order_solution,
                MassMatrix):
        self.clsvof.FCTStep(self.nnz,
                            self.nFreeDOF_global[0],
                            self.weighted_lumped_mass_matrix,
                            self.u_dof_old,
                            high_order_solution,
                            low_order_solution,
                            limited_solution,
                            self.rowptr,  # Row indices for Sparsity Pattern
                            self.colind,  # Column indices for Sparsity Pattern
                            MassMatrix)

    def getNormalReconstruction(self,weighted_mass_matrix):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,weighted_mass_matrix)
        self.clsvof.normalReconstruction(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.offset[0],self.stride[0],
            self.nFreeDOF_global[0], #numDOFs
            self.weighted_lumped_mass_matrix,
            self.rhs_qx,
            self.rhs_qy,
            self.rhs_qz,
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            weighted_mass_matrix)

    def getRhsL2Proj(self):
        self.clsvof.calculateRhsL2Proj(
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.he_for_disc_ICs,
            self.u[0].dof,
            self.offset[0],self.stride[0],
            self.nFreeDOF_global[0], #numDOFs
            self.rhs_l2_proj)

    def getLumpedMassMatrix(self):
        self.lumped_mass_matrix = numpy.zeros(self.u[0].dof.shape,'d')
        self.clsvof.calculateLumpedMassMatrix(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.lumped_mass_matrix,
            self.offset[0],self.stride[0])

    def getResidual(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """

        if self.lumped_mass_matrix is None:
            self.getLumpedMassMatrix()
        if self.getResidualBeforeFirstStep and self.hasVelocityFieldAsFunction:
            self.updateVelocityFieldAsFunction()
        if self.u_dof_old is None:
            # Pass initial condition to u_dof_old. Overwritten if spin up step is taken
            self.u_dof_old = numpy.copy(self.u[0].dof)

        r.fill(0.0)
        self.vofDOFs.fill(0.0)

        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)

        #Dirichlet boundary conditions
        self.numericalFlux.setDirichletValues(self.ebqe)
        #flux boundary conditions
        for t,g in list(self.fluxBoundaryConditionsObjectsDict[0].advectiveFluxBoundaryConditionsDict.items()):
            self.ebqe[('advectiveFlux_bc',0)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
            self.ebqe[('advectiveFlux_bc_flag',0)][t[0],t[1]] = 1

        if self.forceStrongConditions:
              for dofN,g in list(self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.items()):
                  self.u[0].dof[dofN] = g(self.dirichletConditionsForceDOF.DOFBoundaryPointDict[dofN],self.timeIntegration.t)

        min_distance = numpy.zeros(1)
        max_distance = numpy.zeros(1)
        mean_distance = numpy.zeros(1)
        volume_domain = numpy.zeros(1)

        # FREEZE INTERFACE #
        if self.preRedistancingStage==1:
            if self.coefficients.freeze_interface_during_preRedistancing==True:
                for gi in range(len(self.u[0].dof)):
                    if self.interface_locator[gi] == 1.0:
                        self.u[0].dof[gi] = self.u_dof_old[gi]
        # END OF FREEZING INTERFACE #
        else:
            self.interface_locator[:]=0
        #
        self.clsvof.calculateResidual(#element
            self.timeIntegration.dt,
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
            self.mesh.nElements_owned,
            self.coefficients.useMetrics,
            #VRANS start
            self.coefficients.q_vos,
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.degree_polynomial,
            self.u[0].dof,
            self.u_dof_old, #For Backward Euler this is un, for SSP this is the lstage
            self.coefficients.q_v,
            self.coefficients.q_v_old,
            self.timeIntegration.m_tmp[0],
            self.q[('u',0)], #level set
            self.q[('grad(u)',0)], #normal
            self.q[('H(u)',0)], #VOF. Heaviside of level set
            self.q[('mH(u)',0)], #porosity*VOF = (1-vos)*VOF
            self.q['dV'],
            self.q['dV_last'],
            self.q[('cfl',0)],
            self.offset[0],self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            #VRANS start
            self.coefficients.ebqe_vos,
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('u',0)], #level set
            self.ebqe[('grad(u)',0)], #normal
            self.ebqe[('H(u)',0)], #VOF. Heaviside of level set
            self.ebqe[('advectiveFlux',0)],
            # FOR NONLINEAR CLSVOF; i.e., MCorr with VOF
            self.coefficients.timeOrder,
            int(self.timeStage),
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactRedist,
            self.coefficients.lambdaFact,
            # normalization factor
            min_distance,
            max_distance,
            mean_distance,
            volume_domain,
            self.norm_factor_lagged,
            self.VelMax,
            # normal reconstruction
            self.projected_qx_tn,
            self.projected_qy_tn,
            self.projected_qz_tn,
            self.projected_qx_tStar,
            self.projected_qy_tStar,
            self.projected_qz_tStar,
            # To compute H at DOFs
            self.nFreeDOF_global[0], #numDOFs
            self.lumped_mass_matrix,
            self.vofDOFs,
            self.preRedistancingStage,
            self.interface_locator,
            self.coefficients.alpha/self.mesh.elementDiametersArray.min())

        # RELATED TO EIKONAL EQUATION #
        if self.preRedistancingStage == 1:
            # FREEZE INTERFACE #
            if (self.coefficients.freeze_interface_during_preRedistancing==True):
                for gi in range(len(self.u[0].dof)):
                    if self.interface_locator[gi] == 1.0:
                        r[gi] = 0
        # END OF FREEZING INTERFACE #
        else: # RELATED CLSVOF MODEL #
            # Quantities to compute normalization factor
            from proteus.Comm import globalSum, globalMax
            self.min_distance = -globalMax(-min_distance[0])
            self.max_distance = globalMax(max_distance[0])
            self.mean_distance = globalSum(mean_distance[0])
            self.volume_domain = globalSum(volume_domain[0])
            self.mean_distance /= self.volume_domain

            if self.forceStrongConditions:#
                for dofN,g in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.iteritems():
                    r[dofN] = 0

            if self.displayCFL:
                cell_based_cflMax=globalMax(self.q[('cfl',0)].max())*self.timeIntegration.dt
                logEvent("...   Maximum Cell Based CFL = " + str(cell_based_cflMax),level=2)

            if self.stabilization:
                self.stabilization.accumulateSubgridMassHistory(self.q)
            logEvent("Global residual",level=9,data=r)

            self.nonlinear_function_evaluations += 1
            if self.globalResidualDummy is None:
                self.globalResidualDummy = numpy.zeros(r.shape,'d')

    def getJacobian(self,jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)

        self.clsvof.calculateJacobian(#element
            self.timeIntegration.dt,
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
            self.coefficients.useMetrics,
            #VRANS start
            self.coefficients.q_vos,
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.mesh.nodeDiametersArray,
            self.degree_polynomial,
            self.u[0].dof,
            self.u_dof_old, #For Backward/Forward Euler this is un
            self.coefficients.q_v,
            self.q[('cfl',0)],
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            #VRANS start
            self.coefficients.ebqe_vos,
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.csrColumnOffsets_eb[(0,0)],
            self.coefficients.timeOrder,
            int(self.timeStage),
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactRedist,
            self.coefficients.lambdaFact,
            self.preRedistancingStage,
            self.norm_factor_lagged,
            self.coefficients.alpha/self.mesh.elementDiametersArray.min())

        # RELATED TO EIKONAL EQUATION #
        if self.preRedistancingStage == 1:
            # FREEZING INTERFACE #
            if self.coefficients.freeze_interface_during_preRedistancing==True:
                for gi in range(len(self.u[0].dof)):
                    if self.interface_locator[gi] == 1.0:
                        for i in range(self.rowptr[gi], self.rowptr[gi + 1]):
                            if (self.colind[i] == gi):
                                self.nzval[i] = 1.0
                            else:
                                self.nzval[i] = 0.0
            # END OF FREEZING INTERFACE #
        else: # RELATED TO CLSVOF MODEL #
            #Load the Dirichlet conditions directly into residual
            if self.forceStrongConditions:
                scaling = 1.0
                for dofN in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.keys():
                    global_dofN = dofN
                    for i in range(self.rowptr[global_dofN],self.rowptr[global_dofN+1]):
                        if (self.colind[i] == global_dofN):
                            self.nzval[i] = scaling
                        else:
                            self.nzval[i] = 0.0
        #
        logEvent("Jacobian ",level=10,data=jacobian)
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian

    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
                                                 self.q['x'])
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
                                                                                  self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in list(self.advectiveFluxBoundaryConditionsSetterDict.keys())])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
    def estimate_mt(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass
    def updateAfterMeshMotion(self):
        pass
