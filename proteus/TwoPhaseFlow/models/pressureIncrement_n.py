from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from pressureIncrement_p import *
import pressureIncrement_p as physics

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
ct = physics.ct
myTpFlowProblem = physics.myTpFlowProblem
nd = myTpFlowProblem.nd
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
myparams = mparams.pressureIncrement
pparams = params.physical # physical parameters
meshparams = params.mesh

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
he = meshparams.he
triangleFlag = meshparams.triangleFlag
nnx = meshparams.nnx
nny = meshparams.nny
nnz = meshparams.nnz
triangleOptions = meshparams.triangleOptions
parallelPartitioningType = meshparams.parallelPartitioningType
nLayersOfOverlapForParallel = meshparams.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = meshparams.restrictFineSolutionToAllMeshes

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
stepController=FixedStep

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
femSpaces = {0:FESpace['pBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
numericalFluxType = PresInc.NumericalFlux
conservativeFlux=None

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver = LinearSolvers.KSP_petsc4py

linearSmoother = None
if useSuperlu:
    linearSmoother    = None
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'phi_'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
pressure_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
nl_atol_res = pressure_nl_atol_res
tolFac = 0.0
linTolFac = 0.0
l_atol_res = 0.01*pressure_nl_atol_res
maxNonlinearIts = 50

auxiliaryVariables = myparams.auxiliaryVariables
