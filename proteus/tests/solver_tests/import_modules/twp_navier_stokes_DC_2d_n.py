from proteus import *
from twpDC import *
from proteus.default_n import *
from twp_navier_stokes_DC_2d_p import *



class Fixed_dt_controller(proteus.StepControl.Min_dt_controller):
    ''' Class for setting a fixed timestep, dt, from nOptions for the model
    to allow substepping between time intervals in tnlist '''

    def __init__(self,model,nOptions):
        proteus.StepControl.Min_dt_controller.__init__(self,model,nOptions)
        self.dt = nOptions.DT
        self.dt_model_last = None


    def initialize_dt_model(self,t0,tOut):
        self.saveSolution()
        m = self.model.levelModelList[-1]
        self.dt_model = self.dt
        if self.dt_model_last == None:
            self.dt_model_last = self.dt_model
        self.set_dt_allLevels()
        self.substeps = [self.t_model]
        logEvent("Initializing time step on model %s to dt = %12.5e" % (self.model.name,
                                                                   self.dt_model),
            level=1)

    def choose_dt_model(self):
        self.solverFailures=0
        self.errorFailures=0
        self.saveSolution()
        m = self.model.levelModelList[-1]
        self.dt_model = self.dt
        if self.dt_model_last == None:
            self.dt_model_last = self.dt_model
        self.set_dt_allLevels()
        #self.substeps=[self.t_model]
        self.setSubsteps([self.t_model])

    def updateTimeHistory(self,resetFromDOF=False):
        Min_dt_controller.updateTimeHistory(self,resetFromDOF=resetFromDOF)
        self.dt_model_last = self.dt_model

if useBackwardEuler:
    timeIntegration = BackwardEuler#_cfl
    stepController = Min_dt_controller
    #stepController = HeuristicNL_dt_controller
    #nonlinearIterationsFloor = 2
    #nonlinearIterationsCeil=4
    #nonlinearIterationsFloor = 3
    #nonlinearIterationsCeil=4
    #dtNLgrowFactor  = 1.5
    #dtNLreduceFactor= 0.75

    #timeIntegration = NoIntegration
    #stepController = Newton_controller
else:
    timeOrder=2
    timeIntegration = VBDF
    stepController = Fixed_dt_controller #Min_dt_controller
    DT = dt_fixed
    # timeIntegration = FLCBDF
    # stepController  = FLCBDF_controller
    # systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep

    # rtol_u[1] = 1.0e-3
    # rtol_u[2] = 1.0e-3
    # atol_u[1] = 1.0e-2
    # atol_u[2] = 1.0e-2
    # timeOrder = 1
    # class SSPRKwrap(LinearSSPRKintegration):
    #     def __init__(self,vt):
    #         LinearSSPRKintegration.__init__(self,vt,timeOrder,runCFL)
    #         return
    # timeIntegration = SSPRKwrap
    # stepController=Min_dt_RKcontroller
    # systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
    #timeIntegration = VBDF
    # timeIntegration = FLCBDF
    # stepController = Min_dt_controller
    #stepController = FLCBDF_controller_sys
    # rtol_u[1] = 1.0e-2
    # rtol_u[2] = 1.0e-2
    # #rtol_u[3] = 1.0e-2
    # atol_u[1] = 1.0e-2
    # atol_u[2] = 1.0e-2
    # #atol_u[3] = 1.0e-2

femSpaces = {0:basis,
             1:basis,
             2:basis}
numericalFluxType = RANS2P.NumericalFlux
#subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
subgridError = RANS2P.SubgridError(coefficients,nd,lag=False,hFactor=hFactor)
#shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=False)

massLumping = False

fullNewtonFlag = True
multilevelNonlinearSolver = Newton#NS
levelNonlinearSolver = Newton#NS

nonlinearSmoother = None
linearSmoother = None #SimpleNavierStokes2D

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'rans2p_'
#    linearSmoother = StarILU
    linearSmoother = NavierStokes3D_PCD
    linearSolverConvergenceTest = 'r-true'
    parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
    nLayersOfOverlapForParallel = 0
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

tolFac = 0.0
linTolFac = 0.0001
nl_atol_res = 1.0e-8

maxNonlinearIts = 100
maxLineSearches = 0
