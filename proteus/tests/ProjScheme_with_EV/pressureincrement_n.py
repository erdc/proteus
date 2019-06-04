from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from .pressureincrement_p import *

triangleOptions = triangleOptions

femSpaces = {0:pbasis}

stepController=FixedStep

numericalFluxType = PresInc.NumericalFlux
matrix = LinearAlgebraTools.SparseMatrix

linearSmoother    = None
multilevelLinearSolver = LinearSolvers.LU
levelLinearSolver = LinearSolvers.LU

linear_solver_options_prefix = 'phi_'

#fullNewtonFlag = False
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = 0.01*phi_nl_atol_res
tolFac = 0.0
nl_atol_res = phi_nl_atol_res
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0
periodicDirichletConditions=None

#conservativeFlux = {0:'point-eval'} #'point-eval','pwl-bdm-opt'
conservativeFlux=None
