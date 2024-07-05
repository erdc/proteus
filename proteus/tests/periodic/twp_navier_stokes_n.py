from proteus import StepControl, TimeIntegration, NonlinearSolvers, LinearSolvers
from proteus.default_n import *
from twp_navier_stokes_p import *
from proteus import (MeshTools,
                     Context)

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
nnx = ct.nnx
nny = ct.nny
parallelPeriodic=True
if ct.useHex:
    quad=True
triangleFlag=1
#time stepping
runCFL = ct.runCFL
if ct.timeIntegration == "VBDF":
    timeIntegration = TimeIntegration.VBDF
    timeOrder = 2
else:
    timeIntegration = TimeIntegration.BackwardEuler_cfl
    timeOrder = 1
stepController  = StepControl.Min_dt_controller

# mesh options
nLevels = ct.nLevels
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

femSpaces = {0:ct.basis,
             1:ct.basis,
             2:ct.basis}
if nd == 3:
    femSpaces[3] = ct.basis

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients=coefficients,
                                   nd=nd,
                                   lag=ct.ns_lag_subgridError,
                                   hFactor=ct.hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients=coefficients,
                                       nd=nd,
                                       shockCapturingFactor=ct.ns_shockCapturingFactor,
                                       lag=ct.ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
if nd == 2:
    linearSmoother    = LinearSolvers.SimpleNavierStokes2D
elif nd == 3:
    linearSmoother    = LinearSolvers.SimpleNavierStokes3D

matrix = SparseMatrix

multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver      = LinearSolvers.KSP_petsc4py

linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ct.ns_nl_atol_res
nl_atol_res = ct.ns_nl_atol_res
useEisenstatWalker = False#True
maxNonlinearIts = 50
maxLineSearches = 0
#conservativeFlux = {0:'pwl-bdm-opt'}