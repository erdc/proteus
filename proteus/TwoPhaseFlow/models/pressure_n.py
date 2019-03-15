from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from pressure_p import *
<<<<<<< HEAD
=======
import pressure_p as physics
>>>>>>> TwoPhaseFlow

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
<<<<<<< HEAD
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
he = myTpFlowProblem.he
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain
auxVariables = myTpFlowProblem.auxVariables
=======
ct = physics.ct
myTpFlowProblem = physics.myTpFlowProblem
nd = myTpFlowProblem.nd
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
myparams = mparams.pressure
pparams = params.physical # physical parameters
meshparams = params.mesh
>>>>>>> TwoPhaseFlow

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
<<<<<<< HEAD
triangleFlag = myTpFlowProblem.triangleFlag
nnx = myTpFlowProblem.nnx
nny = myTpFlowProblem.nny
nnz = myTpFlowProblem.nnz
triangleOptions = domain.MeshOptions.triangleOptions
=======
he = meshparams.he
triangleFlag = meshparams.triangleFlag
nnx = meshparams.nnx
nny = meshparams.nny
nnz = meshparams.nnz
triangleOptions = meshparams.triangleOptions
parallelPartitioningType = meshparams.parallelPartitioningType
nLayersOfOverlapForParallel = meshparams.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = meshparams.restrictFineSolutionToAllMeshes
>>>>>>> TwoPhaseFlow

# *************************************** #
# ********** TIME INTEGRATION ********** #
# *************************************** #
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
nonlinearSolverConvergenceTest      = 'r'
levelNonlinearSolverConvergenceTest = 'r'

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
numericalFluxType = NumericalFlux.ConstantAdvection_exterior
conservativeFlux=None

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
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
linear_solver_options_prefix = 'pressure_'
linearSolverConvergenceTest         = 'r-true'
maxLineSearches = 0

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
pressure_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
nl_atol_res = pressure_nl_atol_res
tolFac = 0.0
linTolFac = 0.0
l_atol_res = 0.1*pressure_nl_atol_res

<<<<<<< HEAD
# *********************************** #
# ********** AUX VARIABLES ********** #
# *********************************** #
if auxVariables is not None:
    if 'pressure' in auxVariables:
        auxiliaryVariables=auxVariables['pressure']
=======
auxiliaryVariables = myparams.auxiliaryVariables
>>>>>>> TwoPhaseFlow
