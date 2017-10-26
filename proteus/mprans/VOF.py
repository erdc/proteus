import proteus
from proteus.mprans.cVOF import *

class SubgridError(proteus.SubgridError.SGE_base):
    def __init__(self,coefficients,nd):
        proteus.SubgridError.SGE_base.__init__(self,coefficients,nd,lag=False)
    def initializeElementQuadrature(self,mesh,t,cq):
        pass
    def updateSubgridErrorHistory(self,initializationPhase=False):
        pass
    def calculateSubgridError(self,q):
        pass

class ShockCapturing(proteus.ShockCapturing.ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,nStepsToDelay=None):
        proteus.ShockCapturing.ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
        if self.lag:
            logEvent("VOF.ShockCapturing: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay=1
            self.lag=False
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff=[]
        self.numDiff_last=[]
        for ci in range(self.nc):
            self.numDiff.append(cq[('numDiff',ci,ci)])
            self.numDiff_last.append(cq[('numDiff',ci,ci)])
    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay != None and self.nSteps > self.nStepsToDelay:
            logEvent("VOF.ShockCapturing: switched to lagged shock capturing")
            self.lag = True
            self.numDiff_last=[]
            for ci in range(self.nc):
                self.numDiff_last.append(self.numDiff[ci].copy())
        logEvent("VOF: max numDiff %e" % (globalMax(self.numDiff_last[0].max()),))

class NumericalFlux(proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                                        getAdvectiveFluxBoundaryConditions,
                                                                                        getDiffusiveFluxBoundaryConditions)

class RKEV(proteus.TimeIntegration.SSP):
    from proteus import TimeIntegration
    """
    Wrapper for SSPRK time integration using EV

    ... more to come ...
    """
    def __init__(self, transport, timeOrder=1, runCFL=0.1, integrateInterpolationPoints=False):
        TimeIntegration.SSP.__init__(self, transport,integrateInterpolationPoints=integrateInterpolationPoints)
        self.runCFL=runCFL
        self.dtLast=None
        self.dtRatioMax=2.0
        self.isAdaptive=True
        # About the cfl 
        assert transport.coefficients.STABILIZATION_TYPE>0,"SSP method just works for edge based EV methods; i.e., STABILIZATION_TYPE>0"
        assert hasattr(transport,'edge_based_cfl'), "No edge based cfl defined"
        self.cfl = transport.edge_based_cfl
        # Stuff particular for SSP
        self.timeOrder = timeOrder  #order of approximation
        self.nStages = timeOrder  #number of stages total
        self.lstage = 0  #last stage completed
        # storage vectors
        self.u_dof_last = {}
        # per component stage values, list with array at each stage
        self.u_dof_stage = {}
        for ci in range(self.nc):
             if transport.q.has_key(('m',ci)):
                self.u_dof_last[ci] = transport.u[ci].dof.copy()
                self.u_dof_stage[ci] = []
                for k in range(self.nStages+1):                    
                    self.u_dof_stage[ci].append(transport.u[ci].dof.copy())
        
    #def set_dt(self, DTSET):
    #    self.dt = DTSET #  don't update t
    def choose_dt(self):
        maxCFL = 1.0e-6
        maxCFL = max(maxCFL,globalMax(self.cfl.max()))
        self.dt = self.runCFL/maxCFL
        if self.dtLast is None:
            self.dtLast = self.dt
        self.t = self.tLast + self.dt
        if self.dt/self.dtLast  > self.dtRatioMax:
            self.dt = self.dtLast*self.dtRatioMax
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
        assert self.timeOrder in [1,2,3]
        assert self.lstage > 0 and self.lstage <= self.timeOrder
        if self.timeOrder == 3:
            if self.lstage == 1:
                logEvent("First stage of SSP33 method",level=4)
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = self.transport.u[ci].dof
                    # update u_dof_old
                    self.transport.u_dof_old[:] = self.u_dof_stage[ci][self.lstage]
            elif self.lstage == 2:
                logEvent("Second stage of SSP33 method",level=4)
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = self.transport.u[ci].dof
                    self.u_dof_stage[ci][self.lstage] *= 1./4.
                    self.u_dof_stage[ci][self.lstage] += 3./4.*self.u_dof_last[ci]
                    #Update u_dof_old
                    self.transport.u_dof_old[:] = self.u_dof_stage[ci][self.lstage]
            elif self.lstage == 3:
                logEvent("Third stage of SSP33 method",level=4)
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = self.transport.u[ci].dof
                    self.u_dof_stage[ci][self.lstage] *= 2.0/3.0
                    self.u_dof_stage[ci][self.lstage] += 1.0/3.0*self.u_dof_last[ci]
                    # update u_dof_old
                    self.transport.u_dof_old[:] = self.u_dof_last[ci]
                    # update solution to u[0].dof
                    self.transport.u[ci].dof[:] = self.u_dof_stage[ci][self.lstage]
        elif self.timeOrder == 2:
            if self.lstage == 1:
                logEvent("First stage of SSP22 method",level=4)
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = self.transport.u[ci].dof
                    # Update u_dof_old
                    self.transport.u_dof_old[:] = self.transport.u[ci].dof
            elif self.lstage == 2:
                logEvent("Second stage of SSP22 method",level=4)
                for ci in range(self.nc):
                    self.u_dof_stage[ci][self.lstage][:] = self.transport.u[ci].dof
                    self.u_dof_stage[ci][self.lstage][:] *= 1./2.
                    self.u_dof_stage[ci][self.lstage][:] += 1./2.*self.u_dof_last[ci]
                    # update u_dof_old
                    self.transport.u_dof_old[:] = self.u_dof_last[ci]
                    # update solution to u[0].dof
                    self.transport.u[ci].dof[:] = self.u_dof_stage[ci][self.lstage]
        else:
            assert self.timeOrder == 1
            for ci in range(self.nc):
                self.u_dof_stage[ci][self.lstage][:] = self.transport.u[ci].dof[:]

    def initializeTimeHistory(self,resetFromDOF=True):
        """
        Push necessary information into time history arrays
        """
        for ci in range(self.nc):
            self.u_dof_last[ci][:] = self.transport.u[ci].dof[:]
            for k in range(self.nStages):
                self.u_dof_stage[ci][k][:] = self.transport.u[ci].dof[:]

    def updateTimeHistory(self,resetFromDOF=False):
        """
        assumes successful step has been taken
        """
        
        self.t = self.tLast + self.dt
        for ci in range(self.nc):
            self.u_dof_last[ci][:] = self.transport.u[ci].dof[:]
            for k in range(self.nStages):
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
        self.u_dof_stage = {}
        for ci in range(self.nc):
             if self.transport.q.has_key(('m',ci)):
                self.u_dof_stage[ci] = []
                for k in range(self.nStages+1):                    
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
    from proteus.ctransportCoefficients import VOFCoefficientsEvaluate
    from proteus.UnstructuredFMMandFSWsolvers import FMMEikonalSolver,FSWEikonalSolver
    from proteus.NonlinearSolvers import EikonalSolver
    from proteus.ctransportCoefficients import VolumeAveragedVOFCoefficientsEvaluate
    from proteus.cfemIntegrals import copyExteriorElementBoundaryValuesFromElementBoundaryValues
    def __init__(self,
                 LS_model=None,
                 V_model=0,
                 RD_model=None,
                 ME_model=1,
                 EikonalSolverFlag=0,
                 checkMass=True,
                 epsFact=0.0,
                 useMetrics=0.0,
                 sc_uref=1.0,
                 sc_beta=1.0,
                 setParamsFunc=None,
                 movingDomain=False,
                 forceStrongConditions=0,
                 # FOR EDGE BASED EV
                 STABILIZATION_TYPE=0,
                 ENTROPY_TYPE=2, #logarithmic
                 LUMPED_MASS_MATRIX=False,
                 FCT=True,
                 # FOR ENTROPY VISCOSITY
                 cE=1.0,
                 uL=0.0, 
                 uR=1.0,
                 # FOR ARTIFICIAL COMPRESSION
                 cK=1.0,
                 # OUTPUT quantDOFs
                 outputQuantDOFs = False,
                 # NONLINEAR VOF
                 epsFactHeaviside=0.0,
                 epsFactDirac=1.0,
                 epsFactDiffusion=10.0,
                 nonlinearVOF = 0):

        self.epsFactHeaviside=epsFactHeaviside
        self.epsFactDirac=epsFactDirac
        self.epsFactDiffusion=epsFactDiffusion
        self.nonlinearVOF = nonlinearVOF
        self.useMetrics = useMetrics
        self.variableNames=['vof']
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
        self.LS_modelIndex=LS_model
        self.sc_uref=sc_uref
        self.sc_beta=sc_beta
        #mwf added
        self.eikonalSolverFlag = EikonalSolverFlag
        if self.eikonalSolverFlag >= 1: #FMM
            assert self.RD_modelIndex < 0, "no redistance with eikonal solver too"
        self.checkMass = checkMass
        #VRANS
        self.q_porosity = None; self.ebq_porosity = None; self.ebqe_porosity = None
        self.porosity_dof = None
        self.setParamsFunc   = setParamsFunc
        self.flowCoefficients=None
        self.movingDomain=movingDomain
        # EDGE BASED (AND ENTROPY) VISCOSITY 
        self.LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX
        self.STABILIZATION_TYPE=STABILIZATION_TYPE
        self.ENTROPY_TYPE=ENTROPY_TYPE
        self.FCT=FCT
        self.uL=uL
        self.uR=uR
        self.cK=cK
        self.forceStrongConditions=forceStrongConditions
        self.cE=cE
        self.outputQuantDOFs=outputQuantDOFs

    def initializeMesh(self,mesh):
        self.eps = self.epsFact*mesh.h
    def attachModels(self,modelList):
        #self
        self.model = modelList[self.modelIndex]
        #redistanced level set
        if self.RD_modelIndex is not None:
            self.rdModel = modelList[self.RD_modelIndex]
        #level set
        if self.LS_modelIndex is not None:
            self.lsModel = modelList[self.LS_modelIndex]
            self.q_phi = modelList[self.LS_modelIndex].q[('u',0)]
            self.ebqe_phi = modelList[self.LS_modelIndex].ebqe[('u',0)]
            if modelList[self.LS_modelIndex].ebq.has_key(('u',0)):
                self.ebq_phi = modelList[self.LS_modelIndex].ebq[('u',0)]
        else:
            self.ebqe_phi = numpy.zeros(self.model.ebqe[('u',0)].shape,'d')#cek hack, we don't need this        
        #flow model
        #print "flow model index------------",self.flowModelIndex,modelList[self.flowModelIndex].q.has_key(('velocity',0))        
        if self.flowModelIndex is not None:
            if modelList[self.flowModelIndex].q.has_key(('velocity',0)):
                self.q_v = modelList[self.flowModelIndex].q[('velocity',0)]
                self.ebqe_v = modelList[self.flowModelIndex].ebqe[('velocity',0)]
                #self.ux_dof = modelList[solf.flowModelIndex].u[1].dof
                #self.uy_dof = modelList[solf.flowModelIndex].u[2].dof
                #self.uz_dof = modelList[solf.flowModelIndex].u[3].dof
            else:
                self.q_v = modelList[self.flowModelIndex].q[('f',0)]
                self.ebqe_v = modelList[self.flowModelIndex].ebqe[('f',0)]
                #self.ux_dof = modelList[solf.flowModelIndex].u[1].dof
                #self.uy_dof = modelList[solf.flowModelIndex].u[2].dof
                #self.uz_dof = modelList[solf.flowModelIndex].u[3].dof
            if modelList[self.flowModelIndex].ebq.has_key(('velocity',0)):
                self.ebq_v = modelList[self.flowModelIndex].ebq[('velocity',0)]
            else:
                if modelList[self.flowModelIndex].ebq.has_key(('f',0)):
                    self.ebq_v = modelList[self.flowModelIndex].ebq[('f',0)]
        #
        if self.eikonalSolverFlag == 2: #FSW
            self.resDummy = numpy.zeros(self.model.u[0].dof.shape,'d')
            eikonalSolverType = self.FSWEikonalSolver
            self.eikonalSolver = self.EikonalSolver(eikonalSolverType,
                                                    self.model,
                                                    relativeTolerance=0.0,absoluteTolerance=1.0e-12,
                                                    frontTolerance=1.0e-4,#default 1.0e-4
                                                    frontInitType='frontIntersection',#'frontIntersection',#or 'magnitudeOnly'
                                                    useLocalPWLreconstruction = False)
        elif self.eikonalSolverFlag == 1: #FMM
            self.resDummy = numpy.zeros(self.model.u[0].dof.shape,'d')
            eikonalSolverType = self.FMMEikonalSolver
            self.eikonalSolver = self.EikonalSolver(eikonalSolverType,
                                                    self.model,
                                                    frontTolerance=1.0e-4,#default 1.0e-4
                                                    frontInitType='frontIntersection',#'frontIntersection',#or 'magnitudeOnly'
                                                    useLocalPWLreconstruction = False)
        # if self.checkMass:
        #     self.m_pre = Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                              self.model.q[('m',0)],
        #                                              self.model.mesh.nElements_owned)
        #     logEvent("Attach Models VOF: Phase  0 mass after VOF step = %12.5e" % (self.m_pre,),level=2)
        #     self.m_post = Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                              self.model.q[('m',0)],
        #                                              self.model.mesh.nElements_owned)
        #     logEvent("Attach Models VOF: Phase  0 mass after VOF step = %12.5e" % (self.m_post,),level=2)
        #     if self.model.ebqe.has_key(('advectiveFlux',0)):
        #         self.fluxIntegral = Norms.fluxDomainBoundaryIntegral(self.model.ebqe['dS'],
        #                                                              self.model.ebqe[('advectiveFlux',0)],
        #                                                              self.model.mesh)
        #         logEvent("Attach Models VOF: Phase  0 mass conservation after VOF step = %12.5e" % (self.m_post - self.m_pre + self.model.timeIntegration.dt*self.fluxIntegral,),level=2)
        #VRANS
        self.flowCoefficients = modelList[self.flowModelIndex].coefficients
        if hasattr(self.flowCoefficients,'q_porosity'):
            self.q_porosity = self.flowCoefficients.q_porosity            
            if self.STABILIZATION_TYPE>0: # edge based stabilization: EV or smoothness based
                assert hasattr(self.flowCoefficients,'porosity_dof'), 'If STABILIZATION_TYPE>0, the flow model must have porosity_dof'
                self.porosity_dof = self.flowCoefficients.porosity_dof
            else:
                self.porosity_dof = numpy.ones(modelList[self.modelIndex].u[0].dof.shape,'d')
        else:
            # If the flow model doesn't have porosity then set q_porosity=1 and porosity_dof=1
            self.q_porosity = numpy.ones(modelList[self.modelIndex].q[('u',0)].shape,'d')
            self.porosity_dof = numpy.ones(modelList[self.modelIndex].u[0].dof.shape,'d')

            if self.setParamsFunc != None:
                self.setParamsFunc(modelList[self.modelIndex].q['x'],self.q_porosity)
            #
        #
        if hasattr(self.flowCoefficients,'ebq_porosity'):
            self.ebq_porosity = self.flowCoefficients.ebq_porosity
        elif modelList[self.modelIndex].ebq.has_key(('u',0)):
            self.ebq_porosity = numpy.ones(modelList[self.modelIndex].ebq[('u',0)].shape,
                                           'd')
            if self.setParamsFunc != None:
                self.setParamsFunc(modelList[self.modelIndex].ebq['x'],self.ebq_porosity)
            #
        #
        if hasattr(self.flowCoefficients,'ebqe_porosity'):
            self.ebqe_porosity = self.flowCoefficients.ebqe_porosity
        else:
            self.ebqe_porosity = numpy.ones(modelList[self.LS_modelIndex].ebqe[('u',0)].shape,
                                            'd')
            if self.setParamsFunc != None:
                self.setParamsFunc(modelList[self.LS_modelIndex].ebqe['x'],self.ebqe_porosity)
            #
        #
    def initializeElementQuadrature(self,t,cq):
        if self.flowModelIndex == None:
            self.q_v = numpy.ones(cq[('f',0)].shape,'d')
        #VRANS
        self.q_porosity = numpy.ones(cq[('u',0)].shape,'d')

    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        if self.flowModelIndex == None:
            self.ebq_v = numpy.ones(cebq[('f',0)].shape,'d')
        #VRANS
        self.ebq_porosity = numpy.ones(cebq[('u',0)].shape,'d')
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        if self.flowModelIndex == None:
            self.ebqe_v = numpy.ones(cebqe[('f',0)].shape,'d')
        #VRANS
        self.ebqe_porosity = numpy.ones(cebqe[('u',0)].shape,'d')
    def preStep(self,t,firstStep=False):
        # SAVE OLD SOLUTION #
	self.model.u_dof_old[:] = self.model.u[0].dof

        # COMPUTE NEW VELOCITY (if given by user) # 
        if self.model.hasVelocityFieldAsFunction:
            self.model.updateVelocityFieldAsFunction()

        if self.checkMass:
            self.m_pre = Norms.scalarDomainIntegral(self.model.q['dV_last'],
                                                    self.model.q[('m',0)],
                                                    self.model.mesh.nElements_owned)
            logEvent("Phase  0 mass before VOF step = %12.5e" % (self.m_pre,),level=2)
        #     self.m_last = Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                              self.model.timeIntegration.m_last[0],
        #                                              self.model.mesh.nElements_owned)
        #     logEvent("Phase  0 mass before VOF (m_last) step = %12.5e" % (self.m_last,),level=2)
        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        if self.nonlinearVOF > 0:            
            #For visualization, send H(phiHat+u) to the DOFs and store it in quantDOFs
            #This is done via a limited L2Projection
            #Update DOFs of NCLS. Since NCLS, MCorr and CLS live on the same space, this is valid:
            #self.lsModel.u[0].dof += self.model.u[0].dof
            #self.lsModel.q[('u',0)] += self.model.q[('u',0)]
            #self.lsModel.ebqe[('u',0)] += self.model.ebqe[('u',0)]
	    #self.lsModel.q[('grad(u)',0)] += self.model.q[('grad(u)',0)]
	    #self.lsModel.ebqe[('grad(u)',0)] += self.model.ebqe[('grad(u)',0)]            
            #self.lsModel.u[0].dof += self.epsFactDiffusion*self.model.u[0].dof #TMP

            self.model.global_mass_error = 0.0
            self.model.global_L2_interface = 0.0
            self.model.global_H1_interface= 0.0
            self.model.global_L2_Hinterface = 0.0
            self.model.global_H1_Hinterface= 0.0
            self.model.global_L2_u = 0.0
            self.model.global_H1_u = 0.0

            self.model.getRhsL2p(self.model.u[0].dof,
                                 self.lsModel.u[0].dof, #phiHat at dofs
                                 self.lsModel.u_dof_old, #phiExact at dofs
                                 self.model.quantDOFs)
            
            self.model.history.write(repr(self.epsFactDiffusion)+","+
                                     repr(self.model.newton_iterations)+","+
                                     repr(abs(self.model.global_mass_error))+","+
                                     repr(math.sqrt(self.model.global_L2_interface))+","+
                                     repr(math.sqrt(self.model.global_H1_interface))+","+
                                     repr(math.sqrt(self.model.global_L2_Hinterface))+","+
                                     repr(math.sqrt(self.model.global_H1_Hinterface))+","+
                                     repr(math.sqrt(self.model.global_L2_u))+","+
                                     repr(math.sqrt(self.model.global_H1_u))+"\n")
            self.model.history.flush()
        
            #self.model.quantDOFs[:] /= self.model.ML[:]
           

            # DO a nodal projection of H(phiHat)
            from proteus.ctransportCoefficients import smoothedHeaviside
            for i in range (self.model.mesh.nodeDiametersArray.size):
                epsHeaviside = 1.5*self.model.mesh.nodeDiametersArray[i]
                self.model.quantDOFs[i] = smoothedHeaviside(epsHeaviside,
                                                            self.lsModel.u[0].dof[i]
                                                            +self.model.u[0].dof[i])
            self.model.quantDOFs2[:] = self.lsModel.u[0].dof[:] + self.model.u[0].dof[:] #phiHat + u
            
            # Compute norms
            #self.model.
            
        self.model.q['dV_last'][:] = self.model.q['dV']
        if self.checkMass:
            self.m_post = Norms.scalarDomainIntegral(self.model.q['dV'],
                                                     self.model.q[('m',0)],
                                                     self.model.mesh.nElements_owned)
            logEvent("Phase  0 mass after VOF step = %12.5e" % (self.m_post,),level=2)
            #self.fluxIntegral = Norms.fluxDomainBoundaryIntegral(self.model.ebqe['dS'],
            #                                                     self.model.ebqe[('advectiveFlux',0)],
            #                                                     self.model.mesh)
            #logEvent("Phase  0 mass flux boundary integral after VOF step = %12.5e" % (self.fluxIntegral,),level=2)
            #logEvent("Phase  0 mass conservation after VOF step = %12.5e" % (self.m_post - self.m_last + self.model.timeIntegration.dt*self.fluxIntegral,),level=2)
            #divergence = Norms.fluxDomainBoundaryIntegralFromVector(self.model.ebqe['dS'],
            #                                                        self.ebqe_v,
            #                                                        self.model.ebqe['n'],
            #                                                        self.model.mesh)
            #logEvent("Divergence = %12.5e" % (divergence,),level=2)
        copyInstructions = {}
        return copyInstructions
    def updateToMovingDomain(self,t,c):
        #in a moving domain simulation the velocity coming in is already for the moving domain
        pass
    def evaluate(self,t,c):
        #mwf debug
        #print "VOFcoeficients eval t=%s " % t
        if c[('f',0)].shape == self.q_v.shape:
            v = self.q_v
            phi = self.q_phi
            porosity  = self.q_porosity
        elif c[('f',0)].shape == self.ebqe_v.shape:
            v = self.ebqe_v
            phi = self.ebqe_phi
            porosity  = self.ebq_porosity
        elif ((self.ebq_v != None and self.ebq_phi != None) and c[('f',0)].shape == self.ebq_v.shape):
            v = self.ebq_v
            phi = self.ebq_phi
            porosity  = self.ebq_porosity
        else:
            v=None
            phi=None
            porosity=None
        if v != None:
            # self.VOFCoefficientsEvaluate(self.eps,
            #                              v,
            #                              phi,
            #                              c[('u',0)],
            #                              c[('m',0)],
            #                              c[('dm',0,0)],
            #                              c[('f',0)],
            #                              c[('df',0,0)])
            self.VolumeAveragedVOFCoefficientsEvaluate(self.eps,
                                                       v,
                                                       phi,
                                                       porosity,
                                                       c[('u',0)],
                                                       c[('m',0)],
                                                       c[('dm',0,0)],
                                                       c[('f',0)],
                                                       c[('df',0,0)])
        # if self.checkMass:
        #     logEvent("Phase  0 mass in eavl = %12.5e" % (Norms.scalarDomainIntegral(self.model.q['dV'],
        #                                                                        self.model.q[('m',0)],
        #                                                                        self.model.mesh.nElements_owned),),level=2)

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
        self.auxiliaryCallCalculateResidual=False
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
        self.u_dof_old = None

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
        self.edge_based_cfl = numpy.zeros(self.u[0].dof.shape)
        #mesh
        self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,3),'d')
        self.q[('u',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('dV_u',0)] = (1.0/self.mesh.nElements_global)*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('grad(u)',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,self.nSpace_global),'d')
        self.q[('m',0)] = self.q[('u',0)]
        self.q[('m_last',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('mt',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['dV'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q['dV_last'] = -1000*numpy.ones((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('m_tmp',0)] = self.q[('u',0)]
        self.q[('cfl',0)] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.q[('numDiff',0,0)] =  numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element),'d')
        self.ebqe[('u',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('grad(u)',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary,self.nSpace_global),'d')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        # mql. Allow the user to provide functions to define the velocity field
        self.hasVelocityFieldAsFunction = False
        if ('velocityFieldAsFunction') in dir (options): 
            self.velocityFieldAsFunction = options.velocityFieldAsFunction
            self.hasVelocityFieldAsFunction = True

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

        # mql. Some ASSERTS to restrict the combination of the methods
        if self.coefficients.STABILIZATION_TYPE>0:
            assert self.timeIntegration.isSSP==True, "If STABILIZATION_TYPE>0, use RKEV timeIntegration within VOF model"
            cond = 'levelNonlinearSolver' in dir(options) and (options.levelNonlinearSolver==ExplicitLumpedMassMatrix or options.levelNonlinearSolver==ExplicitConsistentMassMatrixForVOF)
            assert cond, "If STABILIZATION_TYPE>0, use levelNonlinearSolver=ExplicitLumpedMassMatrix or ExplicitConsistentMassMatrixForVOF"
        if 'levelNonlinearSolver' in dir(options) and options.levelNonlinearSolver==ExplicitLumpedMassMatrix:
            assert self.coefficients.LUMPED_MASS_MATRIX, "If levelNonlinearSolver=ExplicitLumpedMassMatrix, use LUMPED_MASS_MATRIX=True"
        if self.coefficients.LUMPED_MASS_MATRIX==True:
            cond = self.coefficients.STABILIZATION_TYPE==2 
            assert cond, "Use lumped mass matrix just with: STABILIZATION_TYPE=2 (smoothness based stab.)"
            cond = 'levelNonlinearSolver' in dir(options) and options.levelNonlinearSolver==ExplicitLumpedMassMatrix
            assert cond,"Use levelNonlinearSolver=ExplicitLumpedMassMatrix when the mass matrix is lumped"
        if self.coefficients.FCT==True: 
            cond = self.coefficients.STABILIZATION_TYPE>0, "Use FCT just with STABILIZATION_TYPE>0; i.e., edge based stabilization"
        # END OF ASSERTS 

        #cek adding empty data member for low order numerical viscosity structures here for now
        self.ML=None #lumped mass matrix
        self.MassMatrix=None #consistent mass matrix
        self.MassMatrix_sparseFactor=None
        self.Jacobian_sparseFactor=None
        self.cterm_global=None
        self.cterm_transpose_global=None
        # dL_global and dC_global are not the full matrices but just the CSR arrays containing the non zero entries
        self.low_order_solution=None
        self.dt_times_dC_minus_dL=None
        self.min_u_bc=None
        self.max_u_bc=None
        # Aux quantity at DOFs to be filled by optimized code (MQL)
        self.quantDOFs = numpy.zeros(self.u[0].dof.shape,'d')
        self.quantDOFs2 = numpy.zeros(self.u[0].dof.shape,'d')
        # FOR nonlinear VOF; i.e., MCorr with VOF
        self.phiHat_dof = None
        self.phin = None
        self.calculateResidual = None
        self.calculateJacobian = None
        #FOR limited L2Projection
        self.rhs_mass_correction = None
        self.lumped_L2p = None
        self.limited_L2p = None
        self.consistent_L2p = None
        if self.coefficients.nonlinearVOF>0:
            self.populate_vofModel_with_limited_L2p=False
        else:
            self.populate_vofModel_with_limited_L2p=True

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
        self.vof = cVOF_base(self.nSpace_global,
                             self.nQuadraturePoints_element,
                             self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                             self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                             self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                             self.nElementBoundaryQuadraturePoints_elementBoundary,
                             compKernelFlag)

        self.history = open(self.name+"_lagrange_history.txt","w")
        self.history.write('a'+","+                           
                           'massError'+","+
                           'newtonIterations'+","+
                           'L2norm_interface'+","+
                           'H1norm_interface'+","+
                           'L2norm_Hinterface'+","+
                           'H1norm_Hinterface'+","+
                           'L2norm_u'+","+
                           'H1norm_u'+"\n")
        self.newton_iterations = 0.0
        self.global_mass_error = 0.0
        self.global_L2_interface = 0.0
        self.global_H1_interface= 0.0
        self.global_L2_Hinterface = 0.0
        self.global_H1_Hinterface= 0.0
        self.global_L2_u = 0.0
        self.global_H1_u = 0.0

        
        self.forceStrongConditions=False
        if self.forceStrongConditions:
            self.dirichletConditionsForceDOF = DOFBoundaryConditions(self.u[0].femSpace,dofBoundaryConditionsSetterDict[0],weakDirichletConditions=False)


        if self.movingDomain:
            self.MOVING_DOMAIN=1.0
        else:
            self.MOVING_DOMAIN=0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape,'d')

    def FCTStepL2p(self):
        rowptr, colind, MassMatrix = self.MassMatrix.getCSRrepresentation()
        if (self.limited_L2p is None):
            self.limited_L2p = numpy.zeros(self.u[0].dof.shape,'d')

        self.vof.FCTStepL2p(
            self.nnz, #number of non zero entries 
            len(rowptr)-1, #number of DOFs
            self.ML, #Lumped mass matrix
            self.consistent_L2p, # high order projection
            self.lumped_L2p, #low order projection
            self.limited_L2p,
            rowptr, #Row indices for Sparsity Pattern (convenient for DOF loops)
            colind, #Column indices for Sparsity Pattern (convenient for DOF loops)
            MassMatrix)
        
    def FCTStep(self):
        rowptr, colind, MassMatrix = self.MassMatrix.getCSRrepresentation()
        limited_solution = numpy.zeros(self.u[0].dof.shape)

        self.vof.FCTStep(
                         self.nnz, #number of non zero entries 
                         len(rowptr)-1, #number of DOFs
                         self.ML, #Lumped mass matrix
                         self.timeIntegration.u_dof_stage[0][self.timeIntegration.lstage],#soln
                         self.timeIntegration.u, #high order solution 
                         self.low_order_solution,
                         limited_solution,
                         rowptr, #Row indices for Sparsity Pattern (convenient for DOF loops)
                         colind, #Column indices for Sparsity Pattern (convenient for DOF loops)
                         MassMatrix, 
                         self.dt_times_dC_minus_dL,
                         self.min_u_bc,
                         self.max_u_bc, 
                         self.coefficients.LUMPED_MASS_MATRIX)
        self.timeIntegration.u[:] = limited_solution

    #mwf these are getting called by redistancing classes,
    def calculateCoefficients(self):
        pass

    def updateVelocityFieldAsFunction(self):
        X = {0:self.q[('x')][:,:,0],
             1:self.q[('x')][:,:,1],
             2:self.q[('x')][:,:,2]}
        t = self.timeIntegration.t
        self.coefficients.q_v[...,0] = self.velocityFieldAsFunction[0](X,t)
        self.coefficients.q_v[...,1] = self.velocityFieldAsFunction[1](X,t)
        if (self.nSpace_global==3):
            self.coefficients.q_v[...,2] = self.velocityFieldAsFunction[2](X,t)

        # BOUNDARY
        ebqe_X = {0:self.ebqe['x'][:,:,0],
                  1:self.ebqe['x'][:,:,1],
                  2:self.ebqe['x'][:,:,2]}
        self.coefficients.ebqe_v[...,0] = self.velocityFieldAsFunction[0](ebqe_X,t)
        self.coefficients.ebqe_v[...,1] = self.velocityFieldAsFunction[1](ebqe_X,t)
        if (self.nSpace_global==3):
            self.coefficients.ebqe_v[...,2] = self.velocityFieldAsFunction[2](ebqe_X,t)

    ####################################3
    def getRhsL2p(self,
                        u_dof,
                        phiHat_dof,
                        phiExact_dof,
                        rhs):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """            
        rhs.fill(0.0)

        rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
        (self.global_mass_error,
         self.global_L2_interface,
         self.global_H1_interface,
         self.global_L2_Hinterface,
         self.global_H1_Hinterface,
         self.global_L2_u,
         self.global_H1_u) = self.vof.calculateRhsL2p(#element
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
             self.u[0].femSpace.dofMap.l2g,
             self.mesh.elementDiametersArray,
             u_dof, # This is u_lstage due to update stages in RKEV
             phiHat_dof,
             phiExact_dof,
             self.offset[0],self.stride[0],
             rhs)
    ###############################################
    
    def calculateElementResidual(self):
        if self.globalResidualDummy != None:
            self.getResidual(self.u[0].dof,self.globalResidualDummy)

    def getMassMatrix(self):
        assert self.MassMatrix is not None, "Mass matrix in None, run getResidual first"

    #def setQuantDOFs(self):
    #    from proteus.ctransportCoefficients import smoothedHeaviside
    #    for i in range (self.mesh.nodeDiametersArray.size):
    #        epsHeaviside = self.coefficients.epsFactHeaviside*self.mesh.nodeDiametersArray[i]
    #        phi = self.coefficients.lsModel.u[0].dof[i]
    #        self.quantDOFs[i] = smoothedHeaviside(epsHeaviside,phi)
    #        #self.quantDOFs[i] = 2*smoothedHeaviside(epsHeaviside,phi)-1            
        
    def setMassQuadratureEdgeBasedStabilizationMethods(self):
        if self.rhs_mass_correction is None:
            self.rhs_mass_correction = numpy.zeros(self.u[0].dof.shape,'d')
            self.lumped_L2p = numpy.zeros(self.u[0].dof.shape,'d')
            self.consistent_L2p = numpy.zeros(self.u[0].dof.shape,'d')
        else: 
            self.rhs_mass_correction.fill(0.0)

        mass = self.vof.calculateRhsQuadratureMass(#element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u',0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.psi,
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.offset[0],self.stride[0],
            self.nFreeDOF_global[0], #numDOFs
            self.rhs_mass_correction,
            self.lumped_L2p,
            self.ML,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDiffusion,
            self.phiHat_dof)
        self.mass_history.append([self.timeIntegration.dt,mass])

    def getResidual(self,u,r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        
        if self.coefficients.porosity_dof is None: 
            self.coefficients.porosity_dof = numpy.ones(self.u[0].dof.shape,'d')
        if self.u_dof_old is None:
            # Pass initial condition to u_dof_old
            self.u_dof_old = numpy.copy(self.u[0].dof)
        ########################
        ### COMPUTE C MATRIX ###
        ########################
        if self.cterm_global is None:
            #since we only need cterm_global to persist, we can drop the other self.'s
            self.cterm={}
            self.cterm_a={}
            self.cterm_global={}
            self.cterm_transpose={}
            self.cterm_a_transpose={}
            self.cterm_global_transpose={}
            rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
            nnz = nzval.shape[-1] #number of non-zero entries in sparse matrix
            di = self.q[('grad(u)',0)].copy() #direction of derivative
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
            # SHAPE FUNCTIONS
            self.q[('w',0)] = np.zeros((self.mesh.nElements_global,
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
            self.q[('grad(w)',0)] = np.zeros((self.mesh.nElements_global,
                                              self.nQuadraturePoints_element,
                                              self.nDOF_test_element[0],
                                              self.nSpace_global),
                                             'd')
            self.u[0].femSpace.getBasisGradientValues(self.elementQuadraturePoints,
                                                      self.q['inverse(J)'],
                                                      self.q[('grad(w)',0)])
            self.q[('grad(w)*dV_f',0)] = np.zeros((self.mesh.nElements_global,
                                                   self.nQuadraturePoints_element,
                                                   self.nDOF_test_element[0],
                                                   self.nSpace_global),
                                                  'd')
            cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('u',0)],
                                                          self.q['abs(det(J))'],
                                                          self.q[('grad(w)',0)],
                                                          self.q[('grad(w)*dV_f',0)])
            ##########################
            ### LUMPED MASS MATRIX ###
            ##########################
            #assume a linear mass term
            dm = np.ones(self.q[('u',0)].shape,'d')
            elementMassMatrix = np.zeros((self.mesh.nElements_global,
                                          self.nDOF_test_element[0],
                                          self.nDOF_trial_element[0]),'d')
            cfemIntegrals.updateMassJacobian_weak_lowmem(dm,
                                                         self.q[('w',0)],
                                                         self.q[('w*dV_m',0)],
                                                         elementMassMatrix)
            self.MC_a = nzval.copy()
            self.MassMatrix = SparseMat(self.nFreeDOF_global[0],
                                        self.nFreeDOF_global[0],
                                        nnz,
                                        self.MC_a,
                                        colind,
                                        rowptr)
            cfemIntegrals.zeroJacobian_CSR(self.nnz, self.MassMatrix)
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0,0)],
                                                                      self.csrColumnOffsets[(0,0)],
                                                                      elementMassMatrix,
                                                                      self.MassMatrix)
            self.ML = np.zeros((self.nFreeDOF_global[0],),'d')
            for i in range(self.nFreeDOF_global[0]):
                self.ML[i] = self.MC_a[rowptr[i]:rowptr[i+1]].sum()
            np.testing.assert_almost_equal(self.ML.sum(), 
                                           self.mesh.volume, 
                                           err_msg="Trace of lumped mass matrix should be the domain volume",verbose=True)
            for d in range(self.nSpace_global): #spatial dimensions
                #C matrices
                self.cterm[d] = np.zeros((self.mesh.nElements_global,
                                          self.nDOF_test_element[0],
                                          self.nDOF_trial_element[0]),'d')
                self.cterm_a[d] = nzval.copy()
                #self.cterm_a[d] = numpy.zeros(nzval.size)
                self.cterm_global[d] = SparseMat(self.nFreeDOF_global[0],
                                                 self.nFreeDOF_global[0],
                                                 nnz,
                                                 self.cterm_a[d],
                                                 colind,
                                                 rowptr)
                cfemIntegrals.zeroJacobian_CSR(self.nnz, self.cterm_global[d])
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
                                                                          self.csrRowIndeces[(0,0)],
                                                                          self.csrColumnOffsets[(0,0)],
                                                                          self.cterm[d],
                                                                          self.cterm_global[d])
                #C Transpose matrices
                self.cterm_transpose[d] = np.zeros((self.mesh.nElements_global,
                                                    self.nDOF_test_element[0],
                                                    self.nDOF_trial_element[0]),'d')
                self.cterm_a_transpose[d] = nzval.copy()
                self.cterm_global_transpose[d] = SparseMat(self.nFreeDOF_global[0],
                                                           self.nFreeDOF_global[0],
                                                           nnz,
                                                           self.cterm_a_transpose[d],
                                                           colind,
                                                           rowptr)
                cfemIntegrals.zeroJacobian_CSR(self.nnz, self.cterm_global_transpose[d])
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
                                                                          self.csrRowIndeces[(0,0)],
                                                                          self.csrColumnOffsets[(0,0)],
                                                                          self.cterm_transpose[d],
                                                                          self.cterm_global_transpose[d])
                
        rowptr, colind, Cx = self.cterm_global[0].getCSRrepresentation()
        rowptr, colind, Cy = self.cterm_global[1].getCSRrepresentation()
        if (self.nSpace_global==3):
            rowptr, colind, Cz = self.cterm_global[2].getCSRrepresentation()
        else:
            Cz = numpy.zeros(Cx.shape,'d')
        rowptr, colind, CTx = self.cterm_global_transpose[0].getCSRrepresentation()
        rowptr, colind, CTy = self.cterm_global_transpose[1].getCSRrepresentation()
        if (self.nSpace_global==3):
            rowptr, colind, CTz = self.cterm_global_transpose[2].getCSRrepresentation()
        else:
            CTz = numpy.zeros(CTx.shape,'d')

        # This is dummy. I just care about the csr structure of the sparse matrix
        self.dt_times_dC_minus_dL = np.zeros(Cx.shape,'d')
        self.min_u_bc = numpy.zeros(self.u[0].dof.shape,'d')+1E10
        self.max_u_bc = numpy.zeros(self.u[0].dof.shape,'d')-1E10
        self.low_order_solution = numpy.zeros(self.u[0].dof.shape,'d')
        
        #
        #cek end computationa of cterm_global
        #
        #cek showing mquezada an example of using cterm_global sparse matrix
        #calculation y = c*x where x==1
        #direction=0
        #rowptr, colind, c = self.cterm_global[direction].getCSRrepresentation()
        #y = np.zeros((self.nFreeDOF_global[0],),'d')
        #x = np.ones((self.nFreeDOF_global[0],),'d')
        #ij=0
        #for i in range(self.nFreeDOF_global[0]):
        #    for offset in range(rowptr[i],rowptr[i+1]):
        #        j = colind[offset]
        #        y[i] += c[ij]*x[j]
        #        ij+=1
                        
        r.fill(0.0)
        #Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        #cek can put in logic to skip of BC's don't depend on t or u
        #Dirichlet boundary conditions
        #if hasattr(self.numericalFlux,'setDirichletValues'):
        self.numericalFlux.setDirichletValues(self.ebqe)
        #flux boundary conditions
        for t,g in self.fluxBoundaryConditionsObjectsDict[0].advectiveFluxBoundaryConditionsDict.iteritems():
            self.ebqe[('advectiveFlux_bc',0)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
            self.ebqe[('advectiveFlux_bc_flag',0)][t[0],t[1]] = 1

        if self.forceStrongConditions:
              for dofN,g in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.iteritems():
                  self.u[0].dof[dofN] = g(self.dirichletConditionsForceDOF.DOFBoundaryPointDict[dofN],self.timeIntegration.t)

        degree_polynomial = 1
        try:
            degree_polynomial = self.u[0].femSpace.order
        except:
            pass

        if self.calculateResidual is None:
            if (self.coefficients.STABILIZATION_TYPE == 0): #SUPG
                self.calculateResidual = self.vof.calculateResidual
                self.calculateJacobian = self.vof.calculateJacobian 
            else:
                self.calculateResidual = self.vof.calculateResidual_entropy_viscosity
                self.calculateJacobian = self.vof.calculateMassMatrix

        if self.phiHat_dof is None:
            if self.coefficients.nonlinearVOF > 0:
                cond = self.coefficients.LS_modelIndex is not None
                assert cond, "If nonlinearVOF=True, a NCLS model must be attached"
                self.phiHat_dof = self.coefficients.lsModel.u[0].dof
                self.phin_dof = self.coefficients.lsModel.u_dof_old
                if self.coefficients.nonlinearVOF==1:
                    self.calculateResidual = self.vof.calculateResidual_MCorr_with_VOF
                    self.calculateJacobian = self.vof.calculateJacobian_MCorr_with_VOF
                else:
                    assert self.coefficients.nonlinearVOF==2, 'nonlinearVOF must be 0,1 or 2'
                    self.calculateResidual = self.vof.calculateResidual_MCorr_with_VOF2
                    self.calculateJacobian = self.vof.calculateJacobian_MCorr_with_VOF2
            else:
                self.phiHat_dof = numpy.zeros(self.u[0].dof.shape,'d')
                self.phin_dof = numpy.zeros(self.u[0].dof.shape,'d')

        self.calculateResidual(#element
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
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.shockCapturing.lag,
            self.shockCapturing.shockCapturingFactor,
            self.coefficients.sc_uref,
            self.coefficients.sc_beta,
            #VRANS start
            self.coefficients.q_porosity,
            self.coefficients.porosity_dof, #I need this for edge based methods
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            degree_polynomial,
            self.u[0].dof,
            self.u_dof_old, #For Backward Euler this is un, for SSP this is the lstage
            self.coefficients.q_v,
            self.timeIntegration.m_tmp[0],
            self.q[('u',0)],
            self.timeIntegration.beta_bdf[0],
            self.q['dV'],
            self.q['dV_last'],
            self.q[('cfl',0)],
            self.edge_based_cfl,
            self.shockCapturing.numDiff[0],
            self.shockCapturing.numDiff_last[0],
            self.offset[0],self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            #VRANS start
            self.coefficients.ebqe_porosity,
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.coefficients.ebqe_phi,self.coefficients.epsFact,
            self.ebqe[('u',0)],
            self.ebqe[('advectiveFlux',0)],
            #ENTROPY VISCOSITY and ARTIFICIAL COMRPESSION
            self.coefficients.cE,
            self.coefficients.cK,
            #PARAMETERS FOR LOG BASED ENTROPY FUNCTION
            self.coefficients.uL, 
            self.coefficients.uR,
            #PARAMETERS FOR EDGE VISCOSITY
            len(rowptr)-1, #num of DOFs
            len(Cx), #num of non-zero entries in the sparsity pattern           
            rowptr, #Row indices for Sparsity Pattern (convenient for DOF loops)
            colind, #Column indices for Sparsity Pattern (convenient for DOF loops)
            self.csrRowIndeces[(0,0)], #row indices (convenient for element loops)
            self.csrColumnOffsets[(0,0)], #column indices (convenient for element loops)
            self.csrColumnOffsets_eb[(0, 0)], #indices for boundary terms
            # C matrices
            Cx, 
            Cy,
            Cz,
            CTx,
            CTy,
            CTz,
            self.ML,
            # PARAMETERS FOR 1st or 2nd ORDER MPP METHOD
            self.coefficients.LUMPED_MASS_MATRIX,
            self.coefficients.STABILIZATION_TYPE,
            self.coefficients.ENTROPY_TYPE,
            # FLUX CORRECTED TRANSPORT
            self.low_order_solution,
            self.dt_times_dC_minus_dL, 
            self.min_u_bc,
            self.max_u_bc,
            # FOR NONLINEAR VOF; i.e.,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.phin_dof,
            self.phiHat_dof,
            self.quantDOFs)

        
        if self.forceStrongConditions:#
            for dofN,g in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.iteritems():
                r[dofN] = 0

        if (self.auxiliaryCallCalculateResidual==False):
            edge_based_cflMax=globalMax(self.edge_based_cfl.max())*self.timeIntegration.dt
            cell_based_cflMax=globalMax(self.q[('cfl',0)].max())*self.timeIntegration.dt
            logEvent("...   Current dt = " + str(self.timeIntegration.dt),level=4)
            logEvent("...   Maximum Cell Based CFL = " + str(cell_based_cflMax),level=2)
            logEvent("...   Maximum Edge Based CFL = " + str(edge_based_cflMax),level=2)

        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        logEvent("Global residual",level=9,data=r)
        
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy is None:
            self.globalResidualDummy = numpy.zeros(r.shape,'d')

    def getJacobian(self,jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)

        degree_polynomial = 1
        try:
            degree_polynomial = self.u[0].femSpace.order
        except:
            pass
        
        self.calculateJacobian(#element
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
            self.timeIntegration.alpha_bdf,
            self.shockCapturing.lag,
            self.shockCapturing.shockCapturingFactor,
            #VRANS start
            self.coefficients.q_porosity,
            #VRANS end
            self.u[0].femSpace.dofMap.l2g,
            self.mesh.elementDiametersArray,
            degree_polynomial,
            self.u[0].dof,
            self.coefficients.q_v,
            self.timeIntegration.beta_bdf[0],
            self.q[('cfl',0)],
            self.shockCapturing.numDiff_last[0],
            self.csrRowIndeces[(0,0)],self.csrColumnOffsets[(0,0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            #VRANS start
            self.coefficients.ebqe_porosity,
            #VRANS end
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux_bc_flag',0)],
            self.ebqe[('advectiveFlux_bc',0)],
            self.csrColumnOffsets_eb[(0,0)], 
            self.coefficients.LUMPED_MASS_MATRIX,
            self.coefficients.epsFactHeaviside,
            self.coefficients.epsFactDirac,
            self.coefficients.epsFactDiffusion,
            self.coefficients.cK,
            self.coefficients.uL,
            self.coefficients.uR,
            self.phin_dof,
            self.phiHat_dof)

        #Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            scaling = 1.0#probably want to add some scaling to match non-dirichlet diagonals in linear system
            for dofN in self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.keys():
                    global_dofN = dofN
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
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)
    def estimate_mt(self):
        pass
    def calculateSolutionAtQuadrature(self):
        pass
    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass
    def updateAfterMeshMotion(self):
        pass
