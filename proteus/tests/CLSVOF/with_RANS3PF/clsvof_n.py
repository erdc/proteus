from proteus import *
from proteus.default_n import *
try:
    from .clsvof_p import *
except:
    from clsvof_p import *

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = CLSVOFNewton
fullNewtonFlag = True
updateJacobian = True
timeIntegration = BackwardEuler_cfl

if eps_tolerance_clsvof:
    nl_atol_res = 1E-12
else:
    nl_atol_res=clsvof_nl_atol_res
#
l_atol_res = nl_atol_res
tolFac=0.
maxNonlinearIts = 100000
stepController = Min_dt_controller

femSpaces = {0:pbasis}

#numericalFluxType = CLSVOF.NumericalFlux
#numericalFluxType = DoNothing
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

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

linear_solver_options_prefix = 'clsvof_'
linearSolverConvergenceTest = 'r-true'
