from proteus import *
from threep_navier_stokes_sed_p import *
from dambreak import *

if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
    stepController  = Min_dt_cfl_controller
elif timeDiscretization=='flcbdf':
    timeIntegration = FLCBDF
    #stepController = FLCBDF_controller_sys
    stepController  = Min_dt_cfl_controller
    time_tol = 10.0*ns_sed_nl_atol_res
    atol_u = {0:time_tol,1:time_tol}
    rtol_u = {0:time_tol,1:time_tol}
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller

femSpaces = {0:vbasis,
	     1:vbasis}

massLumping       = False
numericalFluxType = None
conservativeFlux  = None

numericalFluxType = RANS3PSed.NumericalFlux
subgridError = RANS3PSed.SubgridError(coefficients,nd,lag=ns_sed_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS3PSed.ShockCapturing(coefficients,nd,ns_sed_shockCapturingFactor,lag=ns_sed_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None

linearSmoother    = None

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

linear_solver_options_prefix = 'rans3p_sed_'
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest             = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ns_sed_nl_atol_res
nl_atol_res = ns_sed_nl_atol_res
useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
#conservativeFlux = {0:'pwl-bdm-opt'}
#auxiliaryVariables=[pointGauges,lineGauges]
