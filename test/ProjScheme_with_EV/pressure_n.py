from proteus import *
from proteus.default_n import *
from .pressure_p import *


triangleOptions = triangleOptions


femSpaces = {0:pbasis}

stepController=FixedStep

#matrix type
numericalFluxType = NumericalFlux.ConstantAdvection_exterior
matrix = LinearAlgebraTools.SparseMatrix

linear_solver_options_prefix = 'pressure_'

if useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py
    parallelPartitioningType = parallelPartitioningType
    nLayersOfOverlapForParallel = nLayersOfOverlapForParallel
    nonlinearSmoother = None
    linearSmoother    = None

#fullNewtonFlag = False
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = 0.1*pressure_nl_atol_res
tolFac = 0.0
nl_atol_res = pressure_nl_atol_res
maxLineSearches = 0
nonlinearSolverConvergenceTest      = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest         = 'r-true'

periodicDirichletConditions=None

conservativeFlux=None
