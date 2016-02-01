from proteus import StepControl, TimeIntegration, NonlinearSolvers, LinearSolvers
from twp_navier_stokes_p import *
from proteus import Context
ct = Context.get()
from floating_bar import *

timeIntegration = TimeIntegration.BackwardEuler_cfl
stepController  = StepControl.Min_dt_controller

triangleOptions = ct.triangleOptions

femSpaces = {0:ct.basis,
             1:ct.basis,
             2:ct.basis,
             3:ct.basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,ct.nd,lag=ct.ns_lag_subgridError,hFactor=ct.hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,ct.nd,ct.ns_shockCapturingFactor,lag=ct.ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = LinearSolvers.SimpleNavierStokes3D

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

linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*ct.ns_nl_atol_res
nl_atol_res = ct.ns_nl_atol_res
useEisenstatWalker = False#True
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = {0:'pwl-bdm-opt'}
auxiliaryVariables=[ct.bar]
