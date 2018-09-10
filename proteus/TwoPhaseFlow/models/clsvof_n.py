from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from clsvof_p import *

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = CLSVOFNewton
fullNewtonFlag = True
updateJacobian = True
timeIntegration = BackwardEuler_cfl

clsvof_nl_atol_res = max(1.0e-8, 0.001 * ct.he ** 2)
eps_tolerance_clsvof = ct.clsvof_parameters['eps_tolerance_clsvof']
if eps_tolerance_clsvof:
    nl_atol_res = 1E-12
else:
    nl_atol_res=clsvof_nl_atol_res
#
l_atol_res = nl_atol_res
tolFac=0.
maxNonlinearIts = 50
stepController = Min_dt_controller

elementQuadrature = ct.FESpace['elementQuadrature']
elementBoundaryQuadrature = ct.FESpace['elementBoundaryQuadrature']
femSpaces = {0:ct.FESpace['basis']}

#numericalFluxType = CLSVOF.NumericalFlux
#numericalFluxType = DoNothing
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

matrix = SparseMatrix
multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py

if ct.opts.useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'clsvof_'
linearSolverConvergenceTest = 'r-true'
