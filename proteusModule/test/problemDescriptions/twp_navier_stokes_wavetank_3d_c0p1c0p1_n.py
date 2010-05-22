from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_wavetank_3d_p import *
from wavetank3d import *

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,wavetank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,wavetank_quad_order)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=ns_lagSubgridError,delayLagSteps=2)

massLumping = False

shockCapturing = ResGradQuad_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 10

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'pwl-bdm'}
