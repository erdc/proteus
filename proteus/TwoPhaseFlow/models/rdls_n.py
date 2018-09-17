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

ct = physics.ct
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions
params = ct.params

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
triangleFlag = ct.triangleFlag if hasattr(ct,'triangleFlag') else None
nnx = ct.nnx if hasattr(ct,'nnx') else None
nny = ct.nny if hasattr(ct,'nny') else None
nnz = ct.nnz if hasattr(ct,'nnz') else None
triangleOptions = ct.triangleOptions if hasattr(ct,'triangleOptions') and ct.triangleOptions != 'q30DenA' else mesh.triangleOptions
if hasattr(ct, 'nLevels'):
    nLevels = ct.nLevels
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeIntegration = TimeIntegration.NoIntegration
stepController = StepControl.Newton_controller

# ******************************************* #
# ********** FINITE ELEMENT SPACES ********** #
# ******************************************* #
elementQuadrature = ct.FESpace['elementQuadrature']
elementBoundaryQuadrature = ct.FESpace['elementBoundaryQuadrature']
femSpaces = {0: ct.FESpace['lsBasis']}

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
                                     shockCapturingFactor=params.rdls['shockCapturingFactor'],
                                     lag=params.rdls['lag_shockCapturing'])

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
nl_atol_res = max(params.minTol, params.rdls['tolFac']*ct.he**2)
linTolFac = 0.001
l_atol_res = 0.001*nl_atol_res
#
useEisenstatWalker = False#True
tolFac = 0.
maxNonlinearIts = 25
maxLineSearches = 0
