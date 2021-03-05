from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from GN_sw_p import *
from proteus.Transport import Comm

# *********************************************** #
# ********** Read from mySWFlowProblem ********** #
# *********************************************** #
# READ FROM CONTEXT #
runCFL = mySWFlowProblem.cfl
FESpace = mySWFlowProblem.FESpace
he = mySWFlowProblem.he
useSuperlu = mySWFlowProblem.useSuperlu
domain = mySWFlowProblem.domain
SSPOrder = mySWFlowProblem.swe_parameters['SSPOrder']
LUMPED_MASS_MATRIX = mySWFlowProblem.swe_parameters['LUMPED_MASS_MATRIX']
auxiliaryVariables = mySWFlowProblem.auxiliaryVariables

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
if domain is not None:
    triangleFlag = mySWFlowProblem.triangleFlag
    nnx = mySWFlowProblem.nnx
    nny = mySWFlowProblem.nny
    nnz = mySWFlowProblem.nnz
    triangleOptions = domain.MeshOptions.triangleOptions

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeIntegration = GN_SW2DCV.RKEV
timeOrder = SSPOrder
nStagesTime = SSPOrder

# ****************************************** #
# ********** TIME STEP CONTROLLER ********** #
# ****************************************** #
stepController = Min_dt_controller

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
femSpaces = {0: FESpace['basis'],
             1: FESpace['basis'],
             2: FESpace['basis'],
             3: FESpace['basis'],
             4: FESpace['basis'],
             5: FESpace['basis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = Newton
fullNewtonFlag = False  #NOTE: False just if the method is explicit
if (LUMPED_MASS_MATRIX == 1):
    levelNonlinearSolver = ExplicitLumpedMassMatrixShallowWaterEquationsSolver
else:
    levelNonlinearSolver = ExplicitConsistentMassMatrixShallowWaterEquationsSolver

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
try_supg_stabilization = False
subgridError = None
shockCapturing = None
numericalFluxType = GN_SW2DCV.NumericalFlux

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
# change solver for parallel runs
comm = Comm.get()
if comm.size() > 1:
    levelLinearSolver = KSP_petsc4py
    multilevelLinearSolver = KSP_petsc4py
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = 1.0e-7
nl_rtol_res = 0.0
l_atol_res = 1.0e-7
l_rtol_res = 0.0
tolFac = 0.0
maxLineSearches = 0
