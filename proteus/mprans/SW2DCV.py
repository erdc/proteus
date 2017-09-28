import proteus
from proteus.mprans.cSW2DCV import *

class SubgridError(proteus.SubgridError.SGE_base):
    def __init__(self,coefficients,nd,lag=False,nStepsToDelay=0,hFactor=1.0):
        proteus.SubgridError.SGE_base.__init__(self,coefficients,nd,lag)
        self.hFactor=hFactor
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
        if self.lag:
            logEvent("SW2D.SubgridError: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay=1
            self.lag=False
    def initializeElementQuadrature(self,mesh,t,cq):
        import copy
        self.cq=cq
        self.v_last = self.cq[('velocity',0)]
    def updateSubgridErrorHistory(self,initializationPhase=False):
        self.nSteps += 1
        if self.lag:
            self.v_last[:] = self.cq[('velocity',0)]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            logEvent("SW2D.SubgridError: switched to lagged subgrid error")
            self.lag = True
            self.v_last = self.cq[('velocity',0)].copy()
    def calculateSubgridError(self,q):
        pass

class NumericalFlux(proteus.NumericalFlux.ShallowWater_2D):
    hasInterior=False
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 h_eps=1.0e-8,
                 tol_u=1.0e-8):
        proteus.NumericalFlux.ShallowWater_2D.__init__(self,vt,getPointwiseBoundaryConditions,
                                                       getAdvectiveFluxBoundaryConditions,
                                                       getDiffusiveFluxBoundaryConditions,
                                                       getPeriodicBoundaryConditions,
                                                       h_eps,
                                                       tol_u)
        self.penalty_constant = 2.0
        self.includeBoundaryAdjoint=True
        self.boundaryAdjoint_sigma=1.0
        self.hasInterior=False

class ShockCapturing(proteus.ShockCapturing.ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=False,nStepsToDelay=3):
        proteus.ShockCapturing.ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
        if self.lag:
            logEvent("SW2DCV.ShockCapturing: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay=1
            self.lag=False
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff={}
        self.numDiff_last={}
        for ci in range(3):
            self.numDiff[ci] = cq[('numDiff',ci,ci)]
            self.numDiff_last[ci] = cq[('numDiff',ci,ci)]
    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(3):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            logEvent("SW2DCV.ShockCapturing: switched to lagged shock capturing")
            self.lag = True
            for ci in range(3):
                self.numDiff_last[ci] = self.numDiff[ci].copy()
        logEvent("SW2DCV: max numDiff_0 %e numDiff_1 %e numDiff_2 %e" % (globalMax(self.numDiff_last[0].max()),
                                                                    globalMax(self.numDiff_last[1].max()),
                                                                    globalMax(self.numDiff_last[2].max())))

class RKEV(proteus.TimeIntegration.SSP33):
    from proteus import TimeIntegration
    """
    Wrapper for SSPRK time integration using EV

    ... more to come ...
    """
    def __init__(self, transport, timeOrder=1, runCFL=0.1, integrateInterpolationPoints=False):
        BackwardEuler.__init__(self,transport,integrateInterpolationPoints=integrateInterpolationPoints)
        self.trasport = transport
        self.runCFL=runCFL
        self.dtLast=None
        self.isAdaptive=True
        # About the cfl 
        assert hasattr(transport,'edge_based_cfl'), "No edge based cfl defined"
        self.edge_based_cfl = transport.edge_based_cfl
        self.cell_based_cfl = transport.q[('cfl',0)]
        # Stuff particular for SSP33
        self.timeOrder = timeOrder  #order of approximation
        self.nStages = timeOrder  #number of stages total
        self.lstage = 0  #last stage completed
        # storage vectors
        # previous time step mass and solution dof per component
        self.m_last = {}
        #temporarily use this to stash previous solution since m_last used
        #in EV transport models for previous solution value
        self.m_last_save = {}
        self.u_dof_last = {}
        # per component stage values, list with array at each stage
        self.m_stage = {}
        self.u_dof_stage = {}
        for ci in range(self.nc):
             if transport.q.has_key(('m',ci)):
                self.m_last[ci] = transport.q[('m',ci)].copy()
                self.m_last_save[ci] = transport.q[('m',ci)].copy()

                self.u_dof_last[ci] = transport.u[ci].dof.copy()
                self.m_stage[ci] = []
                self.u_dof_stage[ci] = []
                for k in range(self.nStages+1):                    
                    self.m_stage[ci].append(transport.q[('m',ci)].copy())
                    self.u_dof_stage[ci].append(transport.u[ci].dof.copy())
        
    #def set_dt(self, DTSET):
    #    self.dt = DTSET #  don't update t
    def choose_dt(self):
        maxCFL = 1.0e-6
        # COMPUTE edge_based_cfl        
        rowptr_cMatrix, colind_cMatrix, Cx = self.transport.cterm_global[0].getCSRrepresentation()
        rowptr_cMatrix, colind_cMatrix, Cy = self.transport.cterm_global[1].getCSRrepresentation()        
        rowptr_cMatrix, colind_cMatrix, CTx = self.transport.cterm_global_transpose[0].getCSRrepresentation()
        rowptr_cMatrix, colind_cMatrix, CTy = self.transport.cterm_global_transpose[1].getCSRrepresentation()
        numDOFsPerEqn = self.transport.u[0].dof.size
        
        adjusted_maxCFL = self.transport.sw2d.calculateEdgeBasedCFL(
            self.transport.coefficients.g, 
            numDOFsPerEqn,
            self.transport.ML,
            self.transport.u[0].dof, 
            self.transport.u[1].dof, 
            self.transport.u[2].dof,
            self.transport.coefficients.b.dof,
            rowptr_cMatrix, 
            colind_cMatrix, 
            self.transport.hEps, 
            self.transport.hReg, 
            Cx, 
            Cy, 
            CTx, 
            CTy, 
            self.transport.dLow,
            self.runCFL,
            self.transport.edge_based_cfl)

        maxCFL = max(maxCFL,max(adjusted_maxCFL, globalMax(self.edge_based_cfl.max())))
        # maxCFL = max(maxCFL,globalMax(self.cell_based_cfl.max()))
        self.dt = self.runCFL/maxCFL            
        if self.dtLast is None:
            self.dtLast = self.dt
        self.t = self.tLast + self.dt
        # mwf debug
        #print "RKEv max cfl component ci dt dtLast {0} {1} {2} {3}".format(maxCFL,ci,self.dt,self.dtLast)
        self.substeps = [self.t for i in range(self.nStages)] #Manuel is ignoring different time step levels for now
    def initialize_dt(self,t0,tOut,q):
        """
        Modify self.dt
        """
        self.tLast=t0
        self.choose_dt()
        self.t = t0+self.dt

    def setCoefficients(self):
        """
        beta are all 1's here
        mwf not used right now
        """
        self.alpha = numpy.zeros((self.nStages, self.nStages),'d')
        self.dcoefs = numpy.zeros((self.nStages),'d')
        
    def updateStage(self):
        """
        Need to switch to use coefficients
        """
        #mwf debug
        #import pdb
        #pdb.set_trace()
        self.lstage += 1
        assert self.timeOrder in [1,3]
        assert self.lstage > 0 and self.lstage <= self.timeOrder
        if self.timeOrder == 3:
            if self.lstage == 1:
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = numpy.copy(self.transport.u[ci].dof) #no need for .copy?
                    self.m_stage[ci][self.lstage][:] = numpy.copy(self.transport.q[('m',ci)])
                    #needs to be updated for non-scalar equations
                    #this as used as last stage value in EV Transport model
                    #mwf TODO, get rid of m_last here
                    self.m_last[ci] = numpy.copy(self.transport.q[('m',ci)])

            elif self.lstage == 2:
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = numpy.copy(self.transport.u[ci].dof) 
                    self.u_dof_stage[ci][self.lstage] *= 1./4.
                    self.u_dof_stage[ci][self.lstage] += 3./4.*self.u_dof_last[ci]
                    self.m_stage[ci][self.lstage][:] = numpy.copy(self.transport.q[('m',ci)])
                    self.m_stage[ci][self.lstage] *= 1./4.
                    #mwf this has to be fixed
                    #previous stage updated m_last to the stage value
                    #either have another temporary here or have the VOF code use m_stage
                    #instead of m_last
                    self.m_stage[ci][self.lstage] += 3./4.*self.m_last_save[ci] 
                    #mwf TODO get rid of this
                    self.m_last[ci] = numpy.copy(self.m_stage[ci][self.lstage])
            elif self.lstage == 3:
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = numpy.copy(self.transport.u[ci].dof)
                    self.u_dof_stage[ci][self.lstage][:] *= 2.0/3.0
                    self.u_dof_stage[ci][self.lstage][:] += 1.0/3.0*self.u_dof_last[ci]
                    #switch  time history back
                    #mwf TODO this needs to be fixed for multipcomponent
                    self.transport.u[ci].dof[:] = numpy.copy(self.u_dof_stage[ci][self.lstage])
                    self.m_last[ci] = numpy.copy(self.m_last_save[ci])
                    tmp_dof_stage = self.u_dof_stage[ci][self.lstage].copy()
                    self.u_dof_stage[ci][self.lstage][:] = self.u_dof_last[ci]
                    #globalResidualDummy = numpy.zeros(self.u.shape,'d')
                    #self.transport.getResidual(tmp_dof_stage,
                    #                           globalResidualDummy)
                    self.u_dof_stage[ci][self.lstage][:] = tmp_dof_stage
                    #self.setUnkowns(self.u)

        else:
            assert self.timeOrder == 1
            for ci in range(self.nc):
                self.m_stage[ci][self.lstage][:]=self.transport.q[('m',ci)][:]
                self.u_dof_stage[ci][self.lstage][:] = self.transport.u[ci].dof[:]
                            
    def initializeTimeHistory(self,resetFromDOF=True):
        """
        Push necessary information into time history arrays
        """
        for ci in range(self.nc):
            self.m_last[ci][:] = self.transport.q[('m',ci)][:]
            self.u_dof_last[ci][:] = self.transport.u[ci].dof[:]
            self.m_last_save[ci][:] = self.transport.q[('m',ci)][:]
            for k in range(self.nStages):
                self.m_stage[ci][k][:] = self.transport.q[('m',ci)][:]
                self.u_dof_stage[ci][k][:] = self.transport.u[ci].dof[:]
    def updateTimeHistory(self,resetFromDOF=False):
        """
        assumes successful step has been taken
        """
        
        self.t = self.tLast + self.dt
        for ci in range(self.nc):
            self.m_last[ci][:] = self.transport.q[('m',ci)][:]
            self.m_last_save[ci][:] = self.transport.q[('m',ci)][:]
            self.u_dof_last[ci][:] = self.transport.u[ci].dof[:]
            for k in range(self.nStages):
                self.m_stage[ci][k][:]=self.transport.q[('m',ci)][:]
                self.u_dof_stage[ci][k][:] = self.transport.u[ci].dof[:]
        self.lstage=0
        self.dtLast = self.dt
        self.tLast = self.t
    def generateSubsteps(self,tList):
        """
        create list of substeps over time values given in tList. These correspond to stages
        """
        self.substeps = []
        tLast = self.tLast
        for t in tList:
            dttmp = t-tLast
            self.substeps.extend([tLast + dttmp for i in range(self.nStages)])
            tLast = t

    def resetOrder(self,order):
        """
        initialize data structures for stage updges
        """
        self.timeOrder = order  #order of approximation
        self.nStages = order  #number of stages total
        self.lstage = 0  #last stage completed
        # storage vectors
        # per component stage values, list with array at each stage
        self.m_stage = {}
        self.u_dof_stage = {}
        for ci in range(self.nc):
             if self.transport.q.has_key(('m',ci)):
                self.m_stage[ci] = []
                self.u_dof_stage[ci] = []
                for k in range(self.nStages+1):                    
                    self.m_stage[ci].append(self.transport.q[('m',ci)].copy())
                    self.u_dof_stage[ci].append(self.transport.u[ci].dof.copy())
        self.substeps = [self.t for i in range(self.nStages)]            
    def setFromOptions(self,nOptions):
        """
        allow classes to set various numerical parameters
        """
        if 'runCFL' in dir(nOptions):
            self.runCFL = nOptions.runCFL
        flags = ['timeOrder']
        for flag in flags:
            if flag in dir(nOptions):
                val = getattr(nOptions,flag)
                setattr(self,flag,val)
                if flag == 'timeOrder':
                    self.resetOrder(self.timeOrder)

class Coefficients(proteus.TransportCoefficients.TC_base):
    """
    The coefficients for the shallow water equations
    """
    def __init__(self,
                 bathymetry,
                 nu=1.004e-6,
                 g=9.8,
                 nd=2,
                 sd=True,
                 movingDomain=False,
                 useRBLES=0.0,
		 useMetrics=0.0,
                 modelIndex=0, 
                 cE=1.0,
                 LUMPED_MASS_MATRIX=1, 
                 mannings=0.):
        self.bathymetry = bathymetry 
        self.useRBLES=useRBLES
        self.useMetrics=useMetrics
        self.sd=sd
        self.nu = nu
        self.g = g
        self.nd=nd
        self.cE=cE
        self.LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX
        self.mannings=mannings
        self.modelIndex=modelIndex
        mass={}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={}
        if nd==2:
            variableNames=['h','u','v']
            mass= {0:{0:'linear'},
                   1:{0:'linear',1:'linear'},
                   2:{0:'linear',2:'linear'}}
            advection = {0:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear'},
                         1:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear'},
                         2:{0:'nonlinear',
                            1:'nonlinear',
                            2:'nonlinear'}}
            diffusion  = {1:{1:{1:'constant'},2:{2:'constant'}},
                          2:{2:{2:'constant'},1:{1:'constant'}}}
            sdInfo  = {(1,1):(numpy.array([0,1,2],dtype='i'),
                             numpy.array([0,1],dtype='i')),
                       (1,2):(numpy.array([0,0,1],dtype='i'),
                              numpy.array([0],dtype='i')),
                       (2,2):(numpy.array([0,1,2],dtype='i'),
                              numpy.array([0,1],dtype='i')),
                       (2,1):(numpy.array([0,1,1],dtype='i'),
                              numpy.array([1],dtype='i'))}
            potential= {1:{1:'u'},
                        2:{2:'u'}}
            reaction = {1:{0:'linear'},
                        2:{0:'linear'}}
            TC_base.__init__(self,
                             3,
                             mass,
                             advection,
                             diffusion,
                             potential,
                             reaction,
                             hamiltonian,
                             variableNames,
                             sparseDiffusionTensors=sdInfo,
                             useSparseDiffusion = sd,
                             movingDomain=movingDomain)
            self.vectorComponents=[1,2]
    def attachModels(self,modelList):
        self.model = modelList[self.modelIndex]
        #pass
    def initializeMesh(self,mesh):
        x = mesh.nodeArray[:,0]
        y = mesh.nodeArray[:,1]
        if self.bathymetry is None:
            self.b.dof = mesh.nodeArray[:,2].copy()
        else:
            self.b.dof = self.bathymetry[0]([x,y])
        #self.b.dof[:] = 0. #TMP
    def initializeElementQuadrature(self,t,cq):
        pass
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        pass
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        pass
    def updateToMovingDomain(self,t,c):
        pass
    def evaluate(self,t,c):
        pass
    def preStep(self,t,firstStep=False):
        self.model.h_dof_old_old = numpy.copy(self.model.h_dof_old)
        self.model.hu_dof_old_old = numpy.copy(self.model.hu_dof_old)
        self.model.hv_dof_old_old = numpy.copy(self.model.hv_dof_old)
        self.model.h_dof_old = numpy.copy(self.model.u[0].dof)
        self.model.hu_dof_old = numpy.copy(self.model.u[1].dof)
        self.model.hv_dof_old = numpy.copy(self.model.u[2].dof)
    def postStep(self,t,firstStep=False):
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
                 stressTraceBoundaryConditionsSetterDictDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='SW2DCV',
                 reuse_trial_and_test_quadrature=True,
                 sd = True,
                 movingDomain=False,
                 bdyNullSpace=False):
        self.bdyNullSpace=bdyNullSpace
        self.inf_norm_hu=[] #To test 1D well balancing
        self.firstCalculateResidualCall=True
        self.secondCallCalculateResidual=0
        self.postProcessing = False#this is a hack to test the effect of post-processing
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
        #cek hack? give coefficients a bathymetriy array
        import copy
        self.coefficients.b = copy.deepcopy(self.u[0])
        self.coefficients.b.dof.fill(0.0)
        #
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
        #simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q={}
        self.ebq={}
        self.ebq_global={}
        self.ebqe={}
        self.phi_ip={}
        # To compute edge_based_cfl from withing choose_dt of RKEV
        self.edge_based_cfl= numpy.zeros(self.u[0].dof.shape)
        self.dLow=None
        #Old DOFs
        #NOTE (Mql): It is important to link h_dof_old by reference with u[0].dof (and so on).
        # This is because  I need the initial condition to be passed to them as well (before calling calculateResidual). 
        # During preStep I change this and copy the values instead of keeping the reference. 
        self.h_dof_old_old = self.u[0].dof
        self.hu_dof_old_old = self.u[1].dof
        self.hv_dof_old_old = self.u[2].dof
        self.h_dof_old = self.u[0].dof
        self.hu_dof_old = self.u[1].dof
        self.hv_dof_old = self.u[2].dof

        #Vector for mass matrix
        self.check_positivity_water_height=True
        #mesh
        self.h_dof_sge = self.u[0].dof.copy()
        self.hu_dof_sge = self.u[1].dof.copy()
        self.hv_dof_sge = self.u[2].dof.copy()
        self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.ebq_global[('totalFlux',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebq_global[('velocityAverage',0)] = numpy.zeros((self.mesh.nElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.q[('dV_u',0)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dV_u',1)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dV_u',2)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('u',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m',0)] = self.q[('u',0)]
        self.q[('m',1)] = self.q[('u',1)]
        self.q[('m',2)] = self.q[('u',2)]
        self.q[('m_last',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_last',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',1)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',2)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('f',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('velocity',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',0,0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',1,1)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',2,2)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('u',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc_flag',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux_bc',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',1,1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('penalty')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',2,2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('velocity',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('velocity',2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.points_elementBoundaryQuadrature= set()
        self.scalars_elementBoundaryQuadrature= set([('u',ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature= set()
        self.tensors_elementBoundaryQuadrature= set()
        #
        #show quadrature
        #
        logEvent("Dumping quadrature shapes for model %s" % self.name,level=9)
        logEvent("Element quadrature array (q)", level=9)
        for (k,v) in self.q.iteritems(): logEvent(str((k,v.shape)),level=9)
        logEvent("Element boundary quadrature (ebq)",level=9) 
        for (k,v) in self.ebq.iteritems(): logEvent(str((k,v.shape)),level=9)
        logEvent("Global element boundary quadrature (ebq_global)",level=9)
        for (k,v) in self.ebq_global.iteritems(): logEvent(str((k,v.shape)),level=9)
        logEvent("Exterior element boundary quadrature (ebqe)",level=9)
        for (k,v) in self.ebqe.iteritems(): logEvent(str((k,v.shape)),level=9)
        logEvent("Interpolation points for nonlinear diffusion potential (phi_ip)",level=9)
        for (k,v) in self.phi_ip.iteritems(): logEvent(str((k,v.shape)),level=9)
        #
        # allocate residual and Jacobian storage
        #
        #
        # allocate residual and Jacobian storage
        #
        self.elementResidual = [numpy.zeros(
                (self.mesh.nElements_global,
                 self.nDOF_test_element[ci]),
                'd')]
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

        #hEps: this is use to regularize the flux and re-define the dry states
        self.hEps=None
        self.hReg=None
        self.ML=None #lumped mass matrix
        self.MC_global=None #consistent mass matrix
        #Global C Matrices (mql)
        self.cterm_global=None
        self.cterm_transpose_global=None
        # For FCT
        self.low_order_hnp1=None
        self.low_order_hunp1=None
        self.low_order_hvnp1=None
        self.dH_minus_dL=None
        self.muH_minus_muL=None
        # NORMALS 
        self.COMPUTE_NORMALS=1
        self.normalx=None
        self.normaly=None
        self.boundaryIndex=None
        self.reflectingBoundaryConditions=False

        if 'reflecting_BCs' in dir(options) and options.reflecting_BCs:
            self.reflectingBoundaryConditions=True
        # Aux quantity at DOFs to be filled by optimized code (MQL)
        self.quantDOFs = None

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
        if self.postProcessing:
            self.q[('v',0)] = self.tmpvt.q[('v',0)]
            self.ebq[('v',0)] = self.tmpvt.ebq[('v',0)]  
            self.ebq[('w',0)] = self.tmpvt.ebq[('w',0)]
            self.ebq['sqrt(det(g))'] = self.tmpvt.ebq['sqrt(det(g))']
            self.ebq['n'] = self.tmpvt.ebq['n']
            self.ebq[('dS_u',0)] = self.tmpvt.ebq[('dS_u',0)]
            self.ebqe['dS'] = self.tmpvt.ebqe['dS']
            self.ebqe['n'] = self.tmpvt.ebqe['n']
            self.ebq_global['n'] = self.tmpvt.ebq_global['n']
            self.ebq_global['x'] = self.tmpvt.ebq_global['x']
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)  
        logEvent(memory("velocity postprocessor","OneLevelTransport"),level=4)
        #helper for writing out data storage
        from proteus import Archiver
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
        #self.numericalFlux.setDirichletValues(self.ebqe)
        if self.movingDomain:
            self.MOVING_DOMAIN=1.0
        else:
            self.MOVING_DOMAIN=0.0
        #cek hack
        self.movingDomain=False
        self.MOVING_DOMAIN=0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape,'d')
        #cek/ido todo replace python loops in modules with optimized code if possible/necessary
        self.forceStrongConditions=True
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(self.u[cj].femSpace,dofBoundaryConditionsSetterDict[cj],weakDirichletConditions=False)

        compKernelFlag = 0
        #if self.coefficients.useConstantH:
        #    self.elementDiameter = self.mesh.elementDiametersArray.copy()
        #    self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
        #else:
        self.elementDiameter = self.mesh.elementDiametersArray
        print self.nSpace_global," nSpace_global"
        self.sw2d = cSW2DCV_base(self.nSpace_global,
                                   self.nQuadraturePoints_element,
                                   self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                   self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                   self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                   self.nElementBoundaryQuadraturePoints_elementBoundary,
                                   compKernelFlag)

        if 'use_SUPG' in dir(options):
            self.calculateResidual = self.sw2d.calculateResidual_SUPG
            self.calculateJacobian = self.sw2d.calculateJacobian_SUPG
        else:
            self.calculateResidual = self.sw2d.calculateResidual_entropy_viscosity
            if (self.coefficients.LUMPED_MASS_MATRIX):
                self.calculateJacobian = self.sw2d.calculateLumpedMassMatrix
            else:
                self.calculateJacobian = self.sw2d.calculateMassMatrix

    def FCTStep(self):
        #NOTE: this function is meant to be called within the solver
        rowptr, colind, MassMatrix = self.MC_global.getCSRrepresentation()
        # Extract hnp1 from global solution u
        index = range(0,len(self.timeIntegration.u))
        hIndex = index[0::3]
        huIndex = index[1::3]
        hvIndex = index[2::3]
        # create limited solution 
        limited_hnp1 = numpy.zeros(self.h_dof_old.shape)
        limited_hunp1 = numpy.zeros(self.h_dof_old.shape)
        limited_hvnp1 = numpy.zeros(self.h_dof_old.shape)
        # Do some type of limitation 
        
        self.sw2d.FCTStep(self.timeIntegration.dt, 
                          self.nnz, #number of non zero entries 
                          len(rowptr)-1, #number of DOFs
                          self.ML, #Lumped mass matrix
                          self.timeIntegration.u_dof_stage[0][self.timeIntegration.lstage], #hn
                          self.timeIntegration.u_dof_stage[1][self.timeIntegration.lstage], #hun
                          self.timeIntegration.u_dof_stage[2][self.timeIntegration.lstage], #hvn
                          self.coefficients.b.dof,
                          self.timeIntegration.u[hIndex], #high order solution 
                          self.timeIntegration.u[huIndex],
                          self.timeIntegration.u[hvIndex],
                          self.low_order_hnp1, # low order solution 
                          self.low_order_hunp1, 
                          self.low_order_hvnp1, 
                          limited_hnp1, 
                          limited_hunp1,
                          limited_hvnp1,
                          rowptr, #Row indices for Sparsity Pattern (convenient for DOF loops)
                          colind, #Column indices for Sparsity Pattern (convenient for DOF loops)
                          MassMatrix,
                          self.dH_minus_dL,
                          self.muH_minus_muL,
                          self.hEps, 
                          self.hReg,
                          self.coefficients.LUMPED_MASS_MATRIX)
                
        # Pass the post processed hnp1 solution to global solution u
        self.timeIntegration.u[hIndex]  = limited_hnp1
        self.timeIntegration.u[huIndex] = limited_hunp1
        self.timeIntegration.u[hvIndex] = limited_hvnp1

    def getResidual(self,u,r):
        """
        Calculate the element residuals and add in to the global residual
        """

        #COMPUTE hEps
        if self.hEps is None: 
            eps=1E-14
            self.hEps = eps*self.u[0].dof.max()
        #COMPUTE C MATRIX 
        if self.cterm_global is None:
            #since we only need cterm_global to persist, we can drop the other self.'s
            self.cterm={}
            self.cterm_a={}
            self.cterm_global={}
            self.cterm_transpose={}
            self.cterm_a_transpose={}
            self.cterm_global_transpose={}
            #Sparsity pattern for Jacobian
            rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
            nnz = nzval.shape[-1] #number of non-zero entries in sparse matrix

            ###########################################
            ##### SPARSITY PATTERN FOR C MATRICES #####
            ###########################################
            #construct nnz_cMatrix, czval_cMatrix, rowptr_cMatrix, colind_cMatrix C matrix
            nnz_cMatrix = nnz/3/3 #This is always true for the SWEs in 2D
            nzval_cMatrix = numpy.zeros(nnz_cMatrix) #This is enough since the values are filled later
            rowptr_cMatrix = numpy.zeros(self.u[0].dof.size+1,'i') #NOTE that in particular rowptr_cMatrix[0]=0
            colind_cMatrix = numpy.zeros(nnz_cMatrix,'i')
            #fill vector rowptr_cMatrix 
            for i in range(1,rowptr_cMatrix.size):
                rowptr_cMatrix[i] = rowptr_cMatrix[i-1]+(rowptr[3*(i-1)+1]-rowptr[3*(i-1)])/3 
                # = rowptr_cMatrix[i-1] + 1/3*(Number of columns of Jacobian's row 3*(i-1)=0, 3, 6, 9, 12, ... 3*(i-1), ..., 3*n-3)            

            #fill vector colind_cMatrix
            i_cMatrix=0 #ith row of cMatrix
            for i in range(rowptr.size-1): #0 to total num of DOFs (i.e. num of rows of jacobian)
                if (i%3 == 0): # Just consider the rows related to the h variable
                    for j,offset in enumerate(range(rowptr[i],rowptr[i+1])):
                        offset_cMatrix = range(rowptr_cMatrix[i_cMatrix],rowptr_cMatrix[i_cMatrix+1])
                        if (j%3 == 0):
                            colind_cMatrix[offset_cMatrix[j/3]] = colind[offset]/3 
                    i_cMatrix+=1
            # END OF SPARSITY PATTERN FOR C MATRICES 
            
            di = numpy.zeros((self.mesh.nElements_global,
                              self.nQuadraturePoints_element,
                              self.nSpace_global),
                             'd') #direction of derivative
            # JACOBIANS (FOR ELEMENT TRANSFORMATION)
            self.q[('J')] = numpy.zeros((self.mesh.nElements_global,
                                         self.nQuadraturePoints_element,
                                         self.nSpace_global,
                                         self.nSpace_global),
                                        'd')
            self.q[('inverse(J)')] = numpy.zeros((self.mesh.nElements_global,
                                                  self.nQuadraturePoints_element,
                                                  self.nSpace_global,
                                                  self.nSpace_global),
                                                 'd')
            self.q[('det(J)')] = numpy.zeros((self.mesh.nElements_global,
                                              self.nQuadraturePoints_element),
                                             'd')
            self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                             self.q['J'],
                                                             self.q['inverse(J)'],
                                                             self.q['det(J)'])
            self.q['abs(det(J))'] = numpy.abs(self.q['det(J)'])
            # SHAPE FUNCTIONS
            self.q[('w',0)] = numpy.zeros((self.mesh.nElements_global,
                                           self.nQuadraturePoints_element,
                                           self.nDOF_test_element[0]),
                                          'd')
            self.q[('w*dV_m',0)] = self.q[('w',0)].copy()
            self.u[0].femSpace.getBasisValues(self.elementQuadraturePoints, self.q[('w',0)])
            cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u',0)],
                                                 self.q['abs(det(J))'],
                                                 self.q[('w',0)],
                                                 self.q[('w*dV_m',0)])
            #### GRADIENT OF TEST FUNCTIONS 
            self.q[('grad(w)',0)] = numpy.zeros((self.mesh.nElements_global,
                                                 self.nQuadraturePoints_element,
                                                 self.nDOF_test_element[0],
                                                 self.nSpace_global),
                                                'd')
            self.u[0].femSpace.getBasisGradientValues(self.elementQuadraturePoints,
                                                      self.q['inverse(J)'],
                                                      self.q[('grad(w)',0)])
            self.q[('grad(w)*dV_f',0)] = numpy.zeros((self.mesh.nElements_global,
                                                      self.nQuadraturePoints_element,
                                                      self.nDOF_test_element[0],
                                                      self.nSpace_global),
                                                     'd')
            cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('u',0)],
                                                          self.q['abs(det(J))'],
                                                          self.q[('grad(w)',0)],
                                                          self.q[('grad(w)*dV_f',0)])
            #
            #lumped mass matrix
            #
            #assume a linear mass term
            dm = np.ones(self.q[('u',0)].shape,'d')
            elementMassMatrix = np.zeros((self.mesh.nElements_global,
                                          self.nDOF_test_element[0],
                                          self.nDOF_trial_element[0]),'d')
            cfemIntegrals.updateMassJacobian_weak_lowmem(dm,
                                                         self.q[('w',0)],
                                                         self.q[('w*dV_m',0)],
                                                         elementMassMatrix)
            self.MC_a = nzval_cMatrix.copy()
            self.MC_global = SparseMat(self.nFreeDOF_global[0],
                                       self.nFreeDOF_global[0],
                                       nnz_cMatrix,
                                       self.MC_a,
                                       colind_cMatrix,
                                       rowptr_cMatrix)
            cfemIntegrals.zeroJacobian_CSR(nnz_cMatrix, self.MC_global)
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0,0)]/3/3,
                                                                      self.csrColumnOffsets[(0,0)]/3,
                                                                      elementMassMatrix,
                                                                      self.MC_global)
            diamD2 = numpy.sum(self.q['abs(det(J))'][:]*self.elementQuadratureWeights[('u',0)])  
            self.ML = np.zeros((self.nFreeDOF_global[0],),'d')
            self.hReg = np.zeros((self.nFreeDOF_global[0],),'d')
            for i in range(self.nFreeDOF_global[0]):
                self.ML[i] = self.MC_a[rowptr_cMatrix[i]:rowptr_cMatrix[i+1]].sum()
                self.hReg[i] = self.ML[i]/diamD2*self.u[0].dof.max()
            #np.testing.assert_almost_equal(self.ML.sum(), self.mesh.volume, err_msg="Trace of lumped mass matrix should be the domain volume",verbose=True)
            #np.testing.assert_almost_equal(self.ML.sum(), diamD2, err_msg="Trace of lumped mass matrix should be the domain volume",verbose=True)

            for d in range(self.nSpace_global): #spatial dimensions
                #C matrices
                self.cterm[d] = numpy.zeros((self.mesh.nElements_global,
                                             self.nDOF_test_element[0],
                                             self.nDOF_trial_element[0]),'d')
                self.cterm_a[d] = nzval_cMatrix.copy()
                self.cterm_global[d] = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                                    self.nFreeDOF_global[0],
                                                                    nnz_cMatrix,
                                                                    self.cterm_a[d],
                                                                    colind_cMatrix,
                                                                    rowptr_cMatrix)
                cfemIntegrals.zeroJacobian_CSR(nnz_cMatrix, self.cterm_global[d])
                di[:] = 0.0
                di[...,d] = 1.0
                cfemIntegrals.updateHamiltonianJacobian_weak_lowmem(di,
                                                                    self.q[('grad(w)*dV_f',0)],
                                                                    self.q[('w',0)],
                                                                    self.cterm[d]) # int[(di*grad(wj))*wi*dV]
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.csrRowIndeces[(0,0)]/3/3,
                                                                          self.csrColumnOffsets[(0,0)]/3,
                                                                          self.cterm[d],
                                                                          self.cterm_global[d])
                #C Transpose matrices
                self.cterm_transpose[d] = numpy.zeros((self.mesh.nElements_global,
                                                       self.nDOF_test_element[0],
                                                       self.nDOF_trial_element[0]),'d')
                self.cterm_a_transpose[d] = nzval_cMatrix.copy()
                self.cterm_global_transpose[d] = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                                              self.nFreeDOF_global[0],
                                                                              nnz_cMatrix,
                                                                              self.cterm_a_transpose[d],
                                                                              colind_cMatrix,
                                                                              rowptr_cMatrix)
                cfemIntegrals.zeroJacobian_CSR(nnz_cMatrix, self.cterm_global_transpose[d])
                di[:] = 0.0
                di[...,d] = -1.0
                cfemIntegrals.updateAdvectionJacobian_weak_lowmem(di,
                                                                  self.q[('w',0)], 
                                                                  self.q[('grad(w)*dV_f',0)], 
                                                                  self.cterm_transpose[d]) # -int[(-di*grad(wi))*wj*dV]
                cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.l2g[0]['nFreeDOF'],
                                                                          self.l2g[0]['freeLocal'],
                                                                          self.csrRowIndeces[(0,0)]/3/3,
                                                                          self.csrColumnOffsets[(0,0)]/3,
                                                                          self.cterm_transpose[d],
                                                                          self.cterm_global_transpose[d])
        if self.boundaryIndex is None and self.normalx is not None: 
            self.boundaryIndex = []
            for i in range(self.normalx.size):
                if self.normalx[i] != 0 or self.normaly[i] != 0:
                    self.boundaryIndex.append(i)
            self.boundaryIndex = np.array(self.boundaryIndex)
        if self.normalx is None:
            self.normalx = numpy.zeros(self.u[0].dof.shape,'d')
            self.normaly = numpy.zeros(self.u[0].dof.shape,'d')

        if self.quantDOFs is None:
            self.quantDOFs = numpy.zeros(self.u[0].dof.shape,'d')

        rowptr_cMatrix, colind_cMatrix, Cx = self.cterm_global[0].getCSRrepresentation()
        rowptr_cMatrix, colind_cMatrix, Cy = self.cterm_global[1].getCSRrepresentation()        
        rowptr_cMatrix, colind_cMatrix, CTx = self.cterm_global_transpose[0].getCSRrepresentation()
        rowptr_cMatrix, colind_cMatrix, CTy = self.cterm_global_transpose[1].getCSRrepresentation()
        # This is dummy. I just care about the csr structure of the sparse matrix
        self.dH_minus_dL = np.zeros(Cx.shape,'d')
        self.muH_minus_muL = np.zeros(Cx.shape,'d')
        self.low_order_hnp1 = numpy.zeros(self.u[0].dof.shape,'d')
        self.low_order_hunp1 = numpy.zeros(self.u[1].dof.shape,'d')
        self.low_order_hvnp1 = numpy.zeros(self.u[2].dof.shape,'d')

        #Allocate space for dLow (for the first stage in the SSP method)        
        if self.dLow is None: 
            self.dLow = numpy.zeros(Cx.shape,'d')            

        numDOFsPerEqn = self.u[0].dof.size #(mql): I am assuming all variables live on the same FE space
        numNonZeroEntries = len(Cx)
        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        #cek todo put in logic to skip if BC's don't depend on t or u
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
        self.Ct_sge = 4.0
        self.Cd_sge = 144.0
        
        if self.reflectingBoundaryConditions and self.boundaryIndex is not None: 
            self.forceStrongConditions=False
            for dummy, index in enumerate(self.boundaryIndex):
                vx = self.u[1].dof[index]
                vy = self.u[2].dof[index]
                vt = vx*self.normaly[index] - vy*self.normalx[index]
                self.u[1].dof[index] = vt*self.normaly[index]
                self.u[2].dof[index] = -vt*self.normalx[index]

        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                    self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN],self.timeIntegration.t)
        #import pdb
        #pdb.set_trace()
            
        #Make sure that the water height is positive (before computing the residual)
        if (self.check_positivity_water_height==True):
            assert self.u[0].dof.min() >= 0, ("Negative water height: ", self.u[0].dof.min())

        if (self.firstCalculateResidualCall):
            self.timeIntegration.u_dof_stage[0][self.timeIntegration.lstage][:] = self.u[0].dof
            self.timeIntegration.u_dof_stage[1][self.timeIntegration.lstage][:] = self.u[1].dof
            self.timeIntegration.u_dof_stage[2][self.timeIntegration.lstage][:] = self.u[2].dof
            self.firstCalculateResidualCall = False

        self.calculateResidual(
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
            self.elementDiameter,#mesh.elementDiametersArray,
            self.mesh.nElements_global,
            self.coefficients.useRBLES,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.coefficients.nu,
            self.coefficients.g,
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.h_dof_old,
            self.hu_dof_old,
            self.hv_dof_old,
            self.h_dof_old_old,
            self.hu_dof_old_old,
            self.hv_dof_old_old,
            self.timeIntegration.u_dof_stage[0][self.timeIntegration.lstage],
            self.timeIntegration.u_dof_stage[1][self.timeIntegration.lstage],
            self.timeIntegration.u_dof_stage[2][self.timeIntegration.lstage],
            self.coefficients.b.dof,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.h_dof_sge,
            self.hu_dof_sge,
            self.hv_dof_sge,
            self.timeIntegration.m_tmp[0],
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.q[('f',0)],
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.q[('numDiff',0,0)],
            self.q[('numDiff',1,1)], 
            self.q[('numDiff',2,2)], 
            self.shockCapturing.numDiff_last[0],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.coefficients.sdInfo[(1,1)][0],
            self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],
            self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(2,2)][0],
            self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],
            self.coefficients.sdInfo[(2,1)][1],
            self.offset[0],
            self.offset[1],
            self.offset[2],
            self.stride[0],
            self.stride[1],
            self.stride[2],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.q[('velocity',0)],
            self.ebqe[('velocity',0)],
            self.ebq_global[('totalFlux',0)],
            self.elementResidual[0],
            Cx, 
            Cy, 
            CTx, 
            CTy,
            numDOFsPerEqn,
            numNonZeroEntries,
            rowptr_cMatrix,
            colind_cMatrix, 
            self.ML, 
            self.timeIntegration.runCFL,
            self.hEps,
            self.hReg,
            self.q[('u',0)], 
            self.q[('u',1)], 
            self.q[('u',2)],
            self.low_order_hnp1,
            self.low_order_hunp1,
            self.low_order_hvnp1,
            self.dH_minus_dL,
            self.muH_minus_muL,
            self.coefficients.cE, 
            self.coefficients.LUMPED_MASS_MATRIX, 
            self.timeIntegration.dt, 
            self.coefficients.mannings,
            self.quantDOFs, 
            self.secondCallCalculateResidual, 
            self.COMPUTE_NORMALS, 
            self.normalx, 
            self.normaly, 
            self.dLow, 
            self.timeIntegration.lstage)

        self.COMPUTE_NORMALS=0
	if self.forceStrongConditions:#
	    for cj in range(len(self.dirichletConditionsForceDOF)):#
		for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():
                    r[self.offset[cj]+self.stride[cj]*dofN] = 0. #g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN],self.timeIntegration.t)

        if (self.secondCallCalculateResidual==0):
            edge_based_cflMax=globalMax(self.edge_based_cfl.max())*self.timeIntegration.dt
            cell_based_cflMax=globalMax(self.q[('cfl',0)].max())*self.timeIntegration.dt
            logEvent("...   Current dt = " + str(self.timeIntegration.dt),level=4)
            logEvent("...   Maximum Cell Based CFL = " + str(cell_based_cflMax),level=2)
            logEvent("...   Maximum Edge Based CFL = " + str(edge_based_cflMax),level=2)

        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        logEvent("Global residual",level=9,data=r)
        #mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1

    def getJacobian(self,jacobian):
	cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
				       jacobian)
        self.calculateJacobian(
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
            self.elementDiameter,#mesh.elementDiametersArray,
            self.mesh.nElements_global,
            self.coefficients.useRBLES,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.coefficients.nu,
            self.coefficients.g,
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.coefficients.b.dof,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.h_dof_sge,
            self.hu_dof_sge,
            self.hv_dof_sge,
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.stabilization.v_last,
            self.q[('cfl',0)],
            self.shockCapturing.numDiff_last[0],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.coefficients.sdInfo[(1,1)][0],
            self.coefficients.sdInfo[(1,1)][1],
            self.coefficients.sdInfo[(1,2)][0],
            self.coefficients.sdInfo[(1,2)][1],
            self.coefficients.sdInfo[(2,2)][0],
            self.coefficients.sdInfo[(2,2)][1],
            self.coefficients.sdInfo[(2,1)][0],
            self.coefficients.sdInfo[(2,1)][1],
            self.csrRowIndeces[(0,0)],
            self.csrColumnOffsets[(0,0)],
            self.csrRowIndeces[(0,1)],
            self.csrColumnOffsets[(0,1)],
            self.csrRowIndeces[(0,2)],
            self.csrColumnOffsets[(0,2)],
            self.csrRowIndeces[(1,0)],
            self.csrColumnOffsets[(1,0)],
            self.csrRowIndeces[(1,1)],
            self.csrColumnOffsets[(1,1)],
            self.csrRowIndeces[(1,2)],
            self.csrColumnOffsets[(1,2)],
            self.csrRowIndeces[(2,0)],self.csrColumnOffsets[(2,0)],
            self.csrRowIndeces[(2,1)],self.csrColumnOffsets[(2,1)],
            self.csrRowIndeces[(2,2)],self.csrColumnOffsets[(2,2)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc_flag',1)],
            self.ebqe[('advectiveFlux_bc_flag',2)],
            self.ebqe[('diffusiveFlux_bc_flag',1,1)],
            self.ebqe[('diffusiveFlux_bc_flag',2,2)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.ebqe[('advectiveFlux_bc',1)],
            self.ebqe[('advectiveFlux_bc',2)],
            self.numericalFlux.ebqe[('u',1)],
            self.ebqe[('diffusiveFlux_bc',1,1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u',2)],
            self.ebqe[('diffusiveFlux_bc',2,2)],
            self.csrColumnOffsets_eb[(0,0)],
            self.csrColumnOffsets_eb[(0,1)],
            self.csrColumnOffsets_eb[(0,2)],
            self.csrColumnOffsets_eb[(1,0)],
            self.csrColumnOffsets_eb[(1,1)],
            self.csrColumnOffsets_eb[(1,2)],
            self.csrColumnOffsets_eb[(2,0)],
            self.csrColumnOffsets_eb[(2,1)],
            self.csrColumnOffsets_eb[(2,2)], 
            self.timeIntegration.dt)

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
        if self.postProcessing:
            self.tmpvt.calculateElementQuadrature()
        self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t,self.q)
        if self.stabilization is not None:
            self.stabilization.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None:
            self.shockCapturing.initializeElementQuadrature(self.mesh,self.timeIntegration.t,self.q)
    def calculateElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on element boundaries.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.tmpvt.calculateElementBoundaryQuadrature()
        pass
    def calculateExteriorElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.tmpvt.calculateExteriorElementBoundaryQuadrature()
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
        self.fluxBoundaryConditionsObjectsDict = dict([(cj,FluxBoundaryConditions(self.mesh,
                                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                  self.ebqe[('x')],
                                                                                  self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
    def estimate_mt(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        self.h_dof_sge[:] = self.u[0].dof
        self.hu_dof_sge[:] = self.u[1].dof
        self.hv_dof_sge[:] = self.u[2].dof
        # if self.postProcessing:
        #     from proteus.cSW2DCV import calculateVelocityAverage as cva
        #     cva(self.mesh.nExteriorElementBoundaries_global,
        #         self.mesh.exteriorElementBoundariesArray,
        #         self.mesh.nInteriorElementBoundaries_global,
        #         self.mesh.interiorElementBoundariesArray,
        #         self.mesh.elementBoundaryElementsArray,
        #         self.mesh.elementBoundaryLocalElementBoundariesArray,
        #         self.u[1].femSpace.dofMap.l2g,
        #         self.u[1].dof,
        #         self.u[2].dof,
        #         self.u[3].dof,
        #         self.ebq[('v',0)], 
        #         self.ebqe[('velocity',0)],
        #         self.ebq_global[('velocityAverage',0)])
        # self.sw2d.calculateVelocityAverage(self.mesh.nExteriorElementBoundaries_global,
        #                                      self.mesh.exteriorElementBoundariesArray,
        #                                      self.mesh.nInteriorElementBoundaries_global,
        #                                      self.mesh.interiorElementBoundariesArray,
        #                                      self.mesh.elementBoundaryElementsArray,
        #                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
        #                                      self.mesh.nodeArray,
        #                                      self.mesh.elementNodesArray,
        #                                      self.u[0].femSpace.elementMaps.psi_trace,
        #                                      self.u[0].femSpace.elementMaps.grad_psi_trace,
        #                                      self.u[0].femSpace.elementMaps.boundaryNormals,
        #                                      self.u[0].femSpace.elementMaps.boundaryJacobians,
        #                                      self.u[1].femSpace.dofMap.l2g,
        #                                      self.u[1].dof,
        #                                      self.u[2].dof,
        #                                      self.u[3].dof,
        #                                      self.u[1].femSpace.psi_trace,
        #                                      self.ebqe[('velocity',0)],
        #                                      self.ebq_global[('velocityAverage',0)])
        OneLevelTransport.calculateAuxiliaryQuantitiesAfterStep(self)

    def getForce(self,cg,forceExtractionFaces,force,moment):
        pass
        # """
        # Calculate the element residuals and add in to the global residual
        # """

        # #Load the unknowns into the finite element dof
        # #self.timeIntegration.calculateCoefs()
        # #self.timeIntegration.calculateU(u)
        # #self.setUnknowns(self.timeIntegration.u)
        # #cek todo put in logic to skip if BC's don't depend on t or u
        # #hack
        # if self.bcsTimeDependent or not self.bcsSet:
        #     self.bcsSet=True
        #     #Dirichlet boundary conditions
        #     self.numericalFlux.setDirichletValues(self.ebqe)
        #     #Flux boundary conditions
        #     for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
        #         for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
        #             if self.coefficients.advection.has_key(ci):
        #                 self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
        #                 self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
        #         for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
        #             for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
        #                 self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
        #                 self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1

        # # Tag boundaries for force/moment extraction
        # #for cj in range(len(self.dirichletConditionsForceDOF)):
        # #     for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():


	
	# print "SW2DCV Force Faces",len(forceExtractionFaces)

	# #force  = numpy.zeros(3,'d') 
	# #moment = numpy.zeros(3,'d')

        # self.Ct_sge = 4.0
        # self.Cd_sge = 144.0
	# self.C_b    = 10.0
 
        # self.sw2d.calculateForce(#element
        #     self.u[0].femSpace.elementMaps.psi,
        #     self.u[0].femSpace.elementMaps.grad_psi,
        #     self.mesh.nodeArray,
        #     self.mesh.elementNodesArray,
        #     self.elementQuadratureWeights[('u',0)],
        #     self.u[0].femSpace.psi,
        #     self.u[0].femSpace.grad_psi,
        #     self.u[0].femSpace.psi,
        #     self.u[0].femSpace.grad_psi,
        #     self.u[1].femSpace.psi,
        #     self.u[1].femSpace.grad_psi,
        #     self.u[1].femSpace.psi,
        #     self.u[1].femSpace.grad_psi,
        #     #element boundary
        #     self.u[0].femSpace.elementMaps.psi_trace,
        #     self.u[0].femSpace.elementMaps.grad_psi_trace,
        #     self.elementBoundaryQuadratureWeights[('u',0)],
        #     self.u[0].femSpace.psi_trace,
        #     self.u[0].femSpace.grad_psi_trace,
        #     self.u[0].femSpace.psi_trace,
        #     self.u[0].femSpace.grad_psi_trace,
        #     self.u[1].femSpace.psi_trace,
        #     self.u[1].femSpace.grad_psi_trace,
        #     self.u[1].femSpace.psi_trace,
        #     self.u[1].femSpace.grad_psi_trace,
        #     self.u[0].femSpace.elementMaps.boundaryNormals,
        #     self.u[0].femSpace.elementMaps.boundaryJacobians,
        #     #physics
        #     self.mesh.elementDiametersArray,
        #     self.stabilization.hFactor,
        #     self.mesh.nElements_global,
        #     self.coefficients.useRBLES,
        #     self.coefficients.useMetrics,
        #     self.timeIntegration.alpha_bdf,
        #     self.coefficients.epsFact_density,
        #     self.coefficients.epsFact,
        #     self.coefficients.sigma,
        #     self.coefficients.rho_0,
        #     self.coefficients.nu_0,
        #     self.coefficients.rho_1,
        #     self.coefficients.nu_1,
        #     self.Ct_sge,
        #     self.Cd_sge,
        #     self.shockCapturing.shockCapturingFactor,
	#     self.C_b,
        #     self.u[0].femSpace.dofMap.l2g,
        #     self.u[1].femSpace.dofMap.l2g,
        #     self.u[0].dof,
        #     self.u[1].dof,
        #     self.u[2].dof,
        #     self.u[3].dof,
        #     self.coefficients.g,
	#      self.q[('cfl',0)],   # ULTRA UGLY HACK self.q[('rho_0')],
        #     self.coefficients.q_phi,
        #     self.coefficients.q_n,
        #     self.coefficients.q_kappa,
        #     self.timeIntegration.m_tmp[1],
        #     self.timeIntegration.m_tmp[2],
        #     self.timeIntegration.m_tmp[3],
        #     self.q[('f',0)],
        #     self.timeIntegration.beta_bdf[1],
        #     self.timeIntegration.beta_bdf[2],
        #     self.timeIntegration.beta_bdf[3],
        #     self.stabilization.v_last,
        #     self.q[('cfl',0)],
        #     self.q[('numDiff',1,1)], 
        #     self.q[('numDiff',2,2)], 
        #     self.q[('numDiff',3,3)],
        #     self.shockCapturing.numDiff_last[1],
        #     self.shockCapturing.numDiff_last[2],
        #     self.shockCapturing.numDiff_last[3],
        #     self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
        #     self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
        #     self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
        #     self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
        #     self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
        #     self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
        #     self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
        #     self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
        #     self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
        #     self.offset[0],self.offset[1],self.offset[2],self.offset[3],
        #     self.stride[0],self.stride[1],self.stride[2],self.stride[3],
        #     cg, force, moment,
        #     self.mesh.nExteriorElementBoundaries_global,
        #     self.mesh.exteriorElementBoundariesArray,
        #     self.mesh.elementBoundaryElementsArray,
        #     self.mesh.elementBoundaryLocalElementBoundariesArray,
	#     forceExtractionFaces,len(forceExtractionFaces),
        #     self.coefficients.ebqe_phi,
        #     self.coefficients.ebqe_n,
        #     self.coefficients.ebqe_kappa,
        #     self.numericalFlux.isDOFBoundary[0],
        #     self.numericalFlux.isDOFBoundary[1],
        #     self.numericalFlux.isDOFBoundary[2],
        #     self.numericalFlux.isDOFBoundary[3],
        #     self.ebqe[('advectiveFlux_bc_flag',0)],
        #     self.ebqe[('advectiveFlux_bc_flag',1)],
        #     self.ebqe[('advectiveFlux_bc_flag',2)],
        #     self.ebqe[('advectiveFlux_bc_flag',3)],
        #     self.ebqe[('diffusiveFlux_bc_flag',1,1)],
        #     self.ebqe[('diffusiveFlux_bc_flag',2,2)],
        #     self.ebqe[('diffusiveFlux_bc_flag',3,3)],
        #     self.numericalFlux.ebqe[('u',0)],
        #     self.ebqe[('advectiveFlux_bc',0)],
        #     self.ebqe[('advectiveFlux_bc',1)],
        #     self.ebqe[('advectiveFlux_bc',2)],
        #     self.ebqe[('advectiveFlux_bc',3)],
        #     self.numericalFlux.ebqe[('u',1)],
        #     self.ebqe[('diffusiveFlux_bc',1,1)],
        #     self.ebqe[('penalty')],
        #     self.numericalFlux.ebqe[('u',2)],
        #     self.ebqe[('diffusiveFlux_bc',2,2)],
        #     self.numericalFlux.ebqe[('u',3)],
        #     self.ebqe[('diffusiveFlux_bc',3,3)],
        #     self.q[('velocity',0)],
        #     self.ebqe[('velocity',0)],
        #     self.ebq_global[('totalFlux',0)],
        #     self.elementResidual[0])

        # #from mpi4py import MPI	
	# #comm = MPI.COMM_WORLD

	# #tmp1 = numpy.zeros(3,'d')
	# #tmp2 = numpy.zeros(3,'d')	         
	# #comm.Allreduce(force,  tmp1, op=MPI.SUM)     
	# #comm.Allreduce(moment, tmp2, op=MPI.SUM) 
        # #force  [:] = tmp1
	# #moment [:] = tmp2

	# from proteus.flcbdfWrappers import globalSum
        # for i in range(3):
	# 	force[i]  = globalSum(force[i]) 
	# 	moment[i] = globalSum(moment[i]) 

        # #simport time
        # #time.sleep(1)
	# ##comm.Barrier()	
        # ##if self.comm.rank() == 0:
	# #print cg
        # #print "Force and moment in sw2d getForce"
        # #print force 
	# #print moment 
	# ##comm.Barrier()
        # #import time
        # #time.sleep(1)
