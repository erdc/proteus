from importlib import reload
from proteus import *
from proteus.default_n import *
from proteus import defaults
defaults.reset_default_n()
try:
    from . import step2d
except:
    import step2d
reload(step2d)    
try:
    from .step2d import *
except:
    from step2d import *
try:
    from . import twp_navier_stokes_step2d_p
except:
    import twp_navier_stokes_step2d_p
try:
    from .twp_navier_stokes_step2d_p import *
except:
    from twp_navier_stokes_step2d_p import *
from proteus import Context
ct = Context.get()

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

timeIntegration = NoIntegration
stepController = Newton_controller

femSpaces = {0:basis,
             1:basis,
             2:basis}

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

massLumping = False

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None
linearSmoother = None

matrix = SparseMatrix

if usePETSc:    
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'rans2p_'
    schur_solver = ct.opts.schur_solver
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
        # Options for PCD
        # Global KSP options
        from petsc4py import PETSc
        OptDB = PETSc.Options()
        OptDB.clear()
        prefix='rans2p_'
        OptDB.setValue(prefix+'ksp_type', 'fgmres')
        OptDB.setValue(prefix+'ksp_gmres_restart', 300)
        OptDB.setValue(prefix+'ksp_gmres_modifiedgramschmidt', 1)
        OptDB.setValue(prefix+'ksp_pc_side','right')
        OptDB.setValue(prefix+'pc_type', 'fieldsplit')
        OptDB.setValue(prefix+'pc_fieldsplit_type', 'schur')
        OptDB.setValue(prefix+'pc_fieldsplit_schur_fact_type', 'upper')
        OptDB.setValue(prefix+'pc_fieldsplit_schur_precondition', 'user')
        # Velocity block options
        OptDB.setValue(prefix+'fieldsplit_velocity_ksp_type', 'gmres')
        OptDB.setValue(prefix+'fieldsplit_velocity_ksp_gmres_modifiedgramschmidt', 1)
        OptDB.setValue(prefix+'fieldsplit_velocity_ksp_atol', 1e-5)
        OptDB.setValue(prefix+'fieldsplit_velocity_ksp_rtol', 1e-5)
        OptDB.setValue(prefix+'fieldsplit_velocity_ksp_pc_side', 'right')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_u_ksp_type', 'preonly')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_u_pc_type', 'hypre')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_u_pc_hypre_type', 'boomeramg')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_u_pc_hypre_boomeramg_coarsen_type', 'HMIS')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_v_ksp_type', 'preonly')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_v_pc_type', 'hypre')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_v_pc_hypre_type', 'boomeramg')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_v_pc_hypre_boomeramg_coarsen_type', 'HMIS')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_w_ksp_type', 'preonly')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_w_pc_type', 'hypre')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_w_pc_hypre_type', 'boomeramg')
        OptDB.setValue(prefix+'fieldsplit_velocity_fieldsplit_w_pc_hypre_boomeramg_coarsen_type', 'HMIS')
        #PCD Schur Complement options
        OptDB.setValue(prefix+'fieldsplit_pressure_ksp_type', 'preonly')
        OptDB.setValue('innerTPPCDsolver_Qp_visc_ksp_type', 'preonly')
        OptDB.setValue('innerTPPCDsolver_Qp_visc_pc_type', 'lu')
        OptDB.setValue('innerTPPCDsolver_Qp_visc_pc_factor_mat_solver_type', 'superlu_dist')
        OptDB.setValue('innerTPPCDsolver_Qp_dens_ksp_type', 'preonly')
        OptDB.setValue('innerTPPCDsolver_Qp_dens_pc_type', 'lu')
        OptDB.setValue('innerTPPCDsolver_Qp_dens_pc_factor_mat_solver_type', 'superlu_dist')
        OptDB.setValue('innerTPPCDsolver_Ap_rho_ksp_type', 'richardson')
        OptDB.setValue('innerTPPCDsolver_Ap_rho_ksp_max_it', 1)
        #OptDB.setValue('innerTPPCDsolver_Ap_rho_ksp_constant_null_space',1)
        OptDB.setValue('innerTPPCDsolver_Ap_rho_pc_type', 'hypre')
        OptDB.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_type', 'boomeramg')
        OptDB.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_strong_threshold', 0.5)
        OptDB.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_interp_type', 'ext+i-cc')
        OptDB.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_coarsen_type', 'HMIS')
        OptDB.setValue('innerTPPCDsolver_Ap_rho_pc_hypre_boomeramg_agg_nl', 2)
#        linearSmootherOptions = (False, True, True, 3) #(density_scaling, numerical_viscosity, lumped, chebyshev_its)
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
linTolFac = 1.0e-2
l_atol_res = 0.01*ns_nl_atol_res
nl_atol_res = ns_nl_atol_res

maxNonlinearIts = 100
maxLineSearches =0

conservativeFlux=None
