from proteus import *
try:
    from . import cavity2d
except:
    import cavity2d

from proteus.default_n import *
from proteus import defaults
defaults.reset_default_n()
try:
    from .twp_navier_stokes_cavity_2d_p import *
except:
    from twp_navier_stokes_cavity_2d_p import *

from proteus import Context
ct = cavity2d.opts
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
else:
    timeIntegration = BackwardEuler
    stepController = SC_base #FixedStep #SC_base

femSpaces = {0:basis,
             1:basis,
             2:basis}
numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

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
    #schur_solver = 'two_phase_PCD' #'two_phase_PCD' #'selfp_petsc'
    schur_solver = ct.schur_solver
    if schur_solver == 'Qp':
        linearSmoother=NavierStokes3D_Qp
    elif schur_solver == 'petsc_ASM':
        linearSmoother = petsc_ASM
    elif schur_solver == 'two_phase_Qp':
        linearSmoother=NavierStokes_TwoPhaseQp
    elif schur_solver == 'two_phase_LSC':
        linearSmoother=NavierStokes_TwoPhaseLSC
    elif schur_solver == 'two_phase_PCD':
        linearSmoother=NavierStokes_TwoPhasePCD
    elif schur_solver == 'LSC':
        linearSmoother=NavierStokes3D_LSC
    elif schur_solver == 'pcd':
        linearSmoother=NavierStokes3D_PCD
    elif schur_solver == 'selfp_proteus':
        linearSmoother = Schur_Sp
    elif schur_solver == 'selfp_petsc':
        linearSmoother = SimpleNavierStokes2D
    elif schur_solver == 'petsc_LU':
        linearSmoother=petsc_LU
    else:
        raise Exception('invalid solver type')
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
maxLineSearches =0
