from proteus import *
from dambreak import *
from vos_p import *

if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
    stepController  = Min_dt_cfl_controller
elif timeDiscretization=='flcbdf':
    timeIntegration = FLCBDF
    #stepController = FLCBDF_controller
    stepController  = Min_dt_cfl_controller
    time_tol = 10.0*vos_nl_atol_res
    atol_u = {0:time_tol}
    rtol_u = {0:time_tol}
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = VOS3P.NumericalFlux
conservativeFlux  = None
subgridError      = VOS3P.SubgridError(coefficients=coefficients,nd=nd)
shockCapturing    = VOS3P.ShockCapturing(coefficients,nd,shockCapturingFactor=vos_shockCapturingFactor,lag=vos_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

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

linear_solver_options_prefix = 'vos_'
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest         = 'r-true'

tolFac      = 0.0
nl_atol_res = vos_nl_atol_res

linTolFac   = 0.0
l_atol_res = 0.1*vos_nl_atol_res

useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0
