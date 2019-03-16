from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
try:
    from .pressureincrement_p import *
except:
    from pressureincrement_p import *


triangleOptions = triangleOptions

femSpaces = {0:pbasis}

stepController=FixedStep

numericalFluxType = PresInc.NumericalFlux
matrix = LinearAlgebraTools.SparseMatrix

if openTop:
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
else:
    linearSmoother    = LinearSolvers.NavierStokesPressureCorrection # pure neumann laplacian solver
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver = LinearSolvers.KSP_petsc4py

linear_solver_options_prefix = 'phi_'

multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = phi_nl_atol_res
tolFac = 0.0
nl_atol_res = phi_nl_atol_res
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0
periodicDirichletConditions=None

#conservativeFlux = {0:'point-eval'} #'point-eval','pwl-bdm-opt'
conservativeFlux=None
