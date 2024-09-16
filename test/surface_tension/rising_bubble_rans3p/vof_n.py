from proteus import *
try:
    from .risingBubble import *
    from .vof_p import *
except:
    from risingBubble import *
    from vof_p import *
    
if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
    stepController  = Min_dt_cfl_controller
elif timeDiscretization=='flcbdf':
    timeIntegration = FLCBDF
    #stepController = FLCBDF_controller
    stepController  = Min_dt_cfl_controller
    time_tol = 10.0*vof_nl_atol_res
    atol_u = {0:time_tol}
    rtol_u = {0:time_tol}
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller

femSpaces = {0:pbasis}

massLumping       = False
numericalFluxType = VOF3P.NumericalFlux
conservativeFlux  = None
subgridError      = VOF3P.SubgridError(coefficients=coefficients,nd=nd)
shockCapturing    = VOF3P.ShockCapturing(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=vof_lag_shockCapturing)

if EXPLICIT_VOF==True:
    fullNewtonFlag = False
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller
else:
    fullNewtonFlag = True

multilevelNonlinearSolver = Newton
levelNonlinearSolver      = TwoStageNewton

nonlinearSmoother = None
linearSmoother    = None

matrix = SparseMatrix

if useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'vof_'
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'r-true'

tolFac      = 0.0
nl_atol_res = vof_nl_atol_res

linTolFac   = 0.0
l_atol_res = 0.1*vof_nl_atol_res

useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0
