from __future__ import absolute_import
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools,
                     NumericalFlux)
from proteus.mprans import RDLS
import rdls_p as physics

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
ct = physics.ct
myTpFlowProblem = physics.myTpFlowProblem
nd = myTpFlowProblem.nd
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
he = myTpFlowProblem.he
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
triangleFlag = myTpFlowProblem.triangleFlag
nnx = myTpFlowProblem.nnx
nny = myTpFlowProblem.nny
nnz = myTpFlowProblem.nnz
triangleOptions = domain.MeshOptions.triangleOptions
parallelPartitioningType = myTpFlowProblem.parallelPartitioningType
nLayersOfOverlapForParallel = myTpFlowProblem.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = myTpFlowProblem.restrictFineSolutionToAllMeshes

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeIntegration = TimeIntegration.NoIntegration
stepController = StepControl.Newton_controller

# ******************************************* #
# ********** FINITE ELEMENT SPACES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
femSpaces = {0: FESpace['lsBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton
fullNewtonFlag = True
nonlinearSmoother = NonlinearSolvers.NLGaussSeidel
#
levelNonlinearSolverConvergenceTest = 'rits'

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
massLumping = False
numericalFluxType = NumericalFlux.DoNothing
conservativeFlux = None
subgridError = RDLS.SubgridError(coefficients=physics.coefficients,
                                 nd=nd)
shockCapturing = RDLS.ShockCapturing(coefficients=physics.coefficients,
                                     nd=nd,
                                     shockCapturingFactor=mparams.rdls['shockCapturingFactor'],
                                     lag=mparams.rdls['lag_shockCapturing'])

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
linearSmoother = None
multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver = LinearSolvers.KSP_petsc4py
if useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'rdls_'
linearSolverConvergenceTest = 'r-true'


maxNonlinearIts = 25
maxLineSearches = 0
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest = 'r-true'


# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = max(mparams.rdls.minTol, mparams.rdls['tolFac']*he**2)
linTolFac = 0.001
l_atol_res = 0.001*nl_atol_res
#
useEisenstatWalker = False#True
tolFac = 0.
maxNonlinearIts = 25
maxLineSearches = 0
