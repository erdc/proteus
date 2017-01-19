from proteus import *
from kappa_p import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_cfl_controller

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = Kappa.NumericalFlux
conservativeFlux  = None
subgridError      = Kappa.SubgridError(coefficients=coefficients,nd=nd)
shockCapturing    = Kappa.ShockCapturing(coefficients,nd,shockCapturingFactor=kappa_shockCapturingFactor,
                                         lag=kappa_lag_shockCapturing)

fullNewtonFlag  = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
linearSmoother    = None
#printNonlinearSolverInfo = True
matrix = SparseMatrix
if not useOldPETSc and not useSuperlu:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
else:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'kappa_'
levelNonlinearSolverConvergenceTest = 'r'#'rits'
linearSolverConvergenceTest         = 'r'#'rits'

tolFac = 0.0
linTolFac =0.0
l_atol_res = 0.001*kappa_nl_atol_res
nl_atol_res = kappa_nl_atol_res
useEisenstatWalker = True

maxNonlinearIts = 50
maxLineSearches = 0

