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
numericalFluxType = NumericalFlux.Diffusion_IIPG_exterior
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
if ct.opts.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'mesh_'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = max(mparams.moveMeshElastic.minTol, mparams.moveMeshElastic['tolFac']*he**2)
linTolFac = 0.001
l_atol_res = 0.001*nl_atol_res
#
# useEisenstatWalker = False#True
tolFac = 0.
maxNonlinearIts = 4#should be linear
maxLineSearches = 0