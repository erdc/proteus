from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_flume_2d_p import *
from flume import *

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
#timeIntegration = FLCBDF
#stepController  = FLCBDF_controller
#Timeintegration = BackwardEuler
#stepController  = FixedStep
#systemStepControllerType = SplitOperator.Sequential_FixedStep
# timeIntegration = NoIntegration
# stepController  = Newton_controller
#rtol_u[1] = 1.0e-2
#rtol_u[2] = 1.0e-2
#atol_u[1] = 1.0e-2
#atol_u[2] = 1.0e-2

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,flume_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,flume_quad_order)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=ns_lagSubgridError,delayLagSteps=2)

massLumping = False

shockCapturing = ResGradQuad_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 20

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 1.0e-8

conservativeFlux = {0:'pwl-bdm'}
