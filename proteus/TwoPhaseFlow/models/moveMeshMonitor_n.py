from __future__ import absolute_import
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools,
                     NumericalFlux)
import moveMeshMonitor_p as physics

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
mparams = params.Models
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
timeIntegration = TimeIntegration.NoIntegration
# stepController  = StepControl.Min_dt_cfl_controller
# runCFL = ct.opts.cfl

# ******************************************* #
# ********** FINITE ELEMENT SPACES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
femSpaces = {0: FESpace['lsBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = NonlinearSolvers.MoveMeshMonitorNewton
levelNonlinearSolver = NonlinearSolvers.MoveMeshMonitorNewton
nonlinearSmoother = NonlinearSolvers.MoveMeshMonitorNewton
# fullNewtonFlag = True
#
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
# massLumping = False
numericalFluxType = NumericalFlux.Diffusion_SIPG_exterior
conservativeFlux = None
subgridError = None
shockCapturing = None

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
linearSmoother = LinearSolvers.NavierStokesPressureCorrection
multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver = LinearSolvers.KSP_petsc4py
if ct.opts.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'mesh2_'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = max(mparams.moveMeshMonitor.minTol, mparams.moveMeshMonitor['tolFac']*ct.he**2)
linTolFac = 0.
l_atol_res = 0.001*nl_atol_res
#
# useEisenstatWalker = False#True
tolFac = 0.
maxNonlinearIts = 1
maxLineSearches = 0
