from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from clsvof_p import *

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
he = myTpFlowProblem.he
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
triangleFlag = myTpFlowProblem.triangleFlag
nnx = myTpFlowProblem.nnx
nny = myTpFlowProblem.nny
nnz = myTpFlowProblem.nnz
triangleOptions = domain.MeshOptions.triangleOptions

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
runCFL=cfl

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
femSpaces = {0: FESpace['lsBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = CLSVOFNewton
fullNewtonFlag = True
updateJacobian = True
nonlinearSmoother = None

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
numericalFluxType = CLSVOF.NumericalFlux

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = SparseMatrix
linearSmoother = None
multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py
if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU
#
linear_solver_options_prefix = 'clsvof_'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
clsvof_nl_atol_res = max(1.0e-8, 0.001 * he ** 2)
eps_tolerance_clsvof = clsvof_parameters['eps_tolerance_clsvof']
if eps_tolerance_clsvof:
    nl_atol_res = 1E-12
else:
    nl_atol_res=clsvof_nl_atol_res
#
l_atol_res = nl_atol_res
tolFac=0.
maxNonlinearIts = 50
