from __future__ import absolute_import
from proteus import (FemTools,
                     Quadrature,
                     TimeIntegration,
                     NumericalFlux,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
import moveMeshElastic_p as physics
from proteus import Context

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
myparams = mparams.moveMeshMonitor
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
# time stepping
runCFL = cfl
timeIntegration = TimeIntegration.NoIntegration


elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
femSpaces = {0: FESpace['velBasis'],
             1: FESpace['velBasis']}
if nd == 3:
    femSpaces[2] = FESpace['velBasis']

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton
nonlinearSmoother = None
fullNewtonFlag = True
#
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'


# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
massLumping = False
numericalFluxType = NumericalFlux.Stress_IIPG_exterior
conservativeFlux = None
subgridError = None
shockCapturing = None

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
linear_solver_options_prefix = 'mesh_'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = max(myparams.minTol, myparams.tolFac*he**2)
linTolFac = 0.001
l_atol_res = 0.001*nl_atol_res
#
# useEisenstatWalker = False#True
tolFac = 0.
maxNonlinearIts = 4#should be linear
maxLineSearches = 0

auxiliaryVariables = myparams.auxiliaryVariables
