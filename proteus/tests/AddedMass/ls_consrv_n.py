from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools,
                     NumericalFlux)
import ls_consrv_p as physics
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

#time stepping
runCFL = ct.runCFL
timeIntegrator  = TimeIntegration.ForwardIntegrator
timeIntegration = TimeIntegration.NoIntegration

#mesh options
nLevels = ct.nLevels
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes
triangleOptions = mesh.triangleOptions



elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

femSpaces = {0: ct.basis}

subgridError      = None
massLumping       = False
numericalFluxType = NumericalFlux.DoNothing
conservativeFlux  = None
shockCapturing    = None

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = None

matrix = LinearAlgebraTools.SparseMatrix

if ct.useOldPETSc:
    multilevelLinearSolver = LinearSolvers.PETSc
    levelLinearSolver      = LinearSolvers.PETSc
else:
    multilevelLinearSolver = LinearSolvers.KSP_petsc4py
    levelLinearSolver      = LinearSolvers.KSP_petsc4py

if ct.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix = 'mcorr_'
linearSolverConvergenceTest  = 'r-true'

tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*ct.mcorr_nl_atol_res
nl_atol_res = ct.mcorr_nl_atol_res
useEisenstatWalker = False#True

maxNonlinearIts = 50
maxLineSearches = 0
