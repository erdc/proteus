from proteus import *
from ls_p import *

timeIntegration = BackwardEuler_cfl
stepController  = Min_dt_controller

femSpaces = {0:basis}

massLumping       = False
conservativeFlux  = None
numericalFluxType = NCLS.NumericalFlux
subgridError      = NCLS.SubgridError(coefficients,nd)
shockCapturing    = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=ls_lag_shockCapturing)

fullNewtonFlag  = True
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

linear_solver_options_prefix = 'ncls_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*ls_nl_atol_res
nl_atol_res = ls_nl_atol_res
useEisenstatWalker = False

maxNonlinearIts = 100 # 50 Farhad
maxLineSearches = 0
