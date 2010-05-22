from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_bubble_2d_p import *
from bubble import *

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[1] = 1.0e-2
    rtol_u[2] = 1.0e-2
    atol_u[1] = 1.0e-2
    atol_u[2] = 1.0e-2
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[1] = 1.0e-2
    rtol_u[2] = 1.0e-2
    atol_u[1] = 1.0e-2#1.0e-3
    atol_u[2] = 1.0e-2#1.0e-3

noPressureStabilization=False
if spaceOrder==1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis}
    hFactor=1.0
#     femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#                  1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#                  2:C0_AffineQuadraticOnSimplexWithNodalBasis}
#     hFactor=0.5
#     noPressureStabilization=True
#     femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#                  1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#                  2:C0_AffineQuadraticOnSimplexWithNodalBasis}
if spaceOrder==2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 1:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 2:C0_AffineQuadraticOnSimplexWithNodalBasis}
    hFactor=0.5
elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)

elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,bubble_quad_order)

subgridError = None

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=noPressureStabilization)

massLumping = False

shockCapturing = None
shockCapturing = ResGradQuad_SC(coefficients,nd,0.5,lag=True)

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 10
maxLineSearches =25

tolFac = 0.0#0.01

atol = 1.0e-6

matrix = SparseMatrix

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

if usePETSc:    
    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU
    
linearSmoother = GaussSeidel

linTolFac = 1.0e-8

#conservativeFlux = {0:'pwl'}
conservativeFlux = {0:'point-eval'}
