from proteus import *
from proteus.default_n import *
try:
    from .cylinder2d import *
    from .twp_navier_stokes_cylinder_2d_p import *
except:
    from cylinder2d import *
    from twp_navier_stokes_cylinder_2d_p import *

systemStepExact=True

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    #stepController = HeuristicNL_dt_controller
    #nonlinearIterationsFloor = 2
    #nonlinearIterationsCeil=4
    #nonlinearIterationsFloor = 3
    #nonlinearIterationsCeil=4
    #dtNLgrowFactor  = 1.5
    #dtNLreduceFactor= 0.75
    
else:
#     timeIntegration = FLCBDF
#     stepController = FLCBDF_controller_sys
    
    timeOrder=2
    timeIntegration = VBDF 
    stepController = Min_dt_controller
    rtol_u[1] = 1.0e-2
    rtol_u[2] = 1.0e-2
    rtol_u[3] = 1.0e-2
    atol_u[1] = 1.0e-2
    atol_u[2] = 1.0e-2
    atol_u[3] = 1.0e-2

femSpaces = {0:basis,
             1:basis,
             2:basis}

numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients,nd,lag=ns_lag_subgridError,hFactor=hFactor)
shockCapturing = RANS2P.ShockCapturing(coefficients,nd,ns_shockCapturingFactor,lag=ns_lag_shockCapturing)

massLumping = False

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None
linearSmoother = SimpleNavierStokes2D

matrix = SparseMatrix

if usePETSc:    
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'rans2p_'
#    linearSmoother = StarILU
    linearSmoother = None
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

tolFac = 1.0e-10
linTolFac = 0.0
nl_atol_res = 1.0e-10

maxNonlinearIts = 100
maxLineSearches =0
conservativeFlux = {0:'pwl-bdm-opt'}
