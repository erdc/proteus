from pyadh import *
from pyadh.default_n import *
from container import *
from twp_navier_stokes_container_3d_p import *

if useBackwardEuler_ls:
    timeIntegration = BackwardEuler_cfl
    #timeIntegration = BackwardEuler
    stepController = Min_dt_controller
    #stepController = HeuristicNL_dt_controller
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[0] = 1.0e-2
    atol_u[0] = 1.0e-2
                     
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = NavierStokesASGS_velocity_pressure_opt(coefficients,nd,lag=False,delayLagSteps=0,hFactor=1.0)

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior #need weak for parallel and global conservation

massLumping = False

shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 25

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 1.0e-6

nl_atol_res = 1.0e-1

matrix = SparseMatrix

multilevelLinearSolver = PETSc
levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'pwl-bdm'}
