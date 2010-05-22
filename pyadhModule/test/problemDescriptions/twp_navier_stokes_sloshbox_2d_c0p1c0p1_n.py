from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_sloshbox_2d_p import *
from sloshbox import *

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    if timeOrder == 2:
        timeIntegration = VBDF
        stepController = Min_dt_cfl_controller

else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[1] = ns_rtol
    rtol_u[2] = ns_rtol
    atol_u[1] = ns_atol
    atol_u[2] = ns_atol

noPressureStabilization=False
if ns_spaceOrder==1:
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
if ns_spaceOrder==2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 1:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 2:C0_AffineQuadraticOnSimplexWithNodalBasis}
    hFactor=0.5
    noPressureStabilization=False

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

subgridError = None

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=noPressureStabilization)

massLumping = False

#shockCapturing = ResGradQuad_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)
shockCapturing = NavierStokes_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 20
maxLineSearches =25

tolFac = 0.0#0.01

nl_atol_res = ns_atol#0.001*he#1.0e-6

matrix = SparseMatrix

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation
numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior #need weak for parallel and global conservation

if usePETSc:    
    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU
    
linearSmoother = GaussSeidel

linTolFac = 1.0e-8

conservativeFlux = {0:'pwl-bdm'}
#conservativeFlux = {0:'pwl'}
#conservativeFlux = {0:'point-eval'}
