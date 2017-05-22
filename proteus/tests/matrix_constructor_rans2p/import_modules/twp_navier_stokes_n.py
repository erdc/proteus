from proteus import *
from twp_navier_stokes_p import *
from dambreak_Colagrossi_coarse import *

if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
    stepController  = Min_dt_cfl_controller
elif timeDiscretization=='flcbdf':
    timeIntegration = FLCBDF
    #stepController = FLCBDF_controller_sys
    stepController  = Min_dt_cfl_controller
    time_tol = 10.0*ns_nl_atol_res
    atol_u = {1:time_tol,2:time_tol}
    rtol_u = {1:time_tol,2:time_tol}
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller

femSpaces = {0:basis,
	     1:basis,
	     2:basis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None

linearSmoother    = SimpleNavierStokes2D

matrix = SparseMatrix

if useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'rans2p_'
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ns_nl_atol_res
nl_atol_res = ns_nl_atol_res
useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = {0:'pwl-bdm-opt'}
#auxiliaryVariables=[pointGauges,lineGauges]
if schur_solver == 'Qp':
    linearSmoother=NavierStokes3D_Qp
elif schur_solver == 'two_phase_Qp':
    linearSmoother=NavierStokes_TwoPhaseQp
elif schur_solver == 'two_phase_LSC':
    linearSmoother=NavierStokes_TwoPhaseLSC
elif schur_solver == 'two_phase_PCD':
    linearSmoother=NavierStokes_TwoPhasePCD
elif schur_solver == 'LSC':
    linearSmoother=NavierStokes3D_LSC
elif schur_solver == 'pcd':
    linearSmoother=NavierStokes3D_PCD
elif schur_solver == 'selfp':
    linearSmoother=SimpleNavierStokes3D
elif schur_solver == 'petsc_LU':
    linearSmoother=petsc_LU
else:
    raise Exception, 'invalid solver type'
