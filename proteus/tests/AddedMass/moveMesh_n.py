from proteus.default_n import *
from proteus import (FemTools,
                     Quadrature,
                     TimeIntegration,
                     NumericalFlux,
                     NonlinearSolvers,
                     LinearSolvers)
import moveMesh_p as physics
from proteus import Context
ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

# time stepping
runCFL = ct.runCFL
timeIntegration = TimeIntegration.NoIntegration

# mesh options
nLevels = ct.nLevels
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes
triangleOptions = mesh.triangleOptions



elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

femSpaces = {0: ct.basis,
             1: ct.basis}
if nd == 3:
    femSpaces[2] = ct.basis

massLumping       = False
numericalFluxType = NumericalFlux.Stress_IIPG_exterior
conservativeFlux  = None

subgridError = None
shockCapturing = None

fullNewtonFlag = True
multilevelNonlinearSolver  = NonlinearSolvers.Newton
levelNonlinearSolver       = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother = None

matrix = SparseMatrix

if ct.useOldPETSc:
    multilevelLinearSolver = LinearSolvers.PETSc
    levelLinearSolver      = LinearSolvers.PETSc
else:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver      = LinearSolvers.KSP_petsc4py

if ct.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix = 'mesh_'
linearSmoother = None
linearSolverConvergenceTest = 'r-true'

tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*ct.mesh_nl_atol_res
nl_atol_res = ct.mesh_nl_atol_res
maxNonlinearIts = 4#should be linear
maxLineSearches = 0

