from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_sloshbox_3d_p import *
from sloshbox3d import *

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    if timeOrder == 2:
        timeIntegration = VBDF
        stepController = Min_dt_cfl_controller

#    timeIntegration = FLCBDF
#    stepController = FLCBDF_controller_sys
#    rtol_u[1] = 1.0e-2
#    rtol_u[2] = 1.0e-2
#    atol_u[1] = 1.0e-2#1.0e-3
#    atol_u[2] = 1.0e-2#1.0e-3
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[1] = 1.0e-2
    rtol_u[2] = 1.0e-2
    rtol_u[3] = 1.0e-2
    atol_u[1] = 1.0e-2
    atol_u[2] = 1.0e-2
    atol_u[3] = 1.0e-2

noPressureStabilization=False
if spaceOrder==1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis,
                 3:C0_AffineLinearOnSimplexWithNodalBasis}
    hFactor=1.0
if spaceOrder==2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 1:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 2:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 3:C0_AffineQuadraticOnSimplexWithNodalBasis}
    hFactor=0.5
elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

subgridError = None

if useOpt:
    subgridError = NavierStokesASGS_velocity_pressure_opt(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=noPressureStabilization)
    shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=True)
else:
    subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=noPressureStabilization)
    shockCapturing = NavierStokes_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)
massLumping = False



numericalFluxType = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 25
maxLineSearches =50

tolFac = 0.0#0.01

nl_atol_res = 1.0e-4#0.001*he

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

conservativeFlux = None
conservativeFlux = {0:'pwl-bdm'}
#conservativeFlux = {0:'pwl'}
#conservativeFlux = {0:'point-eval'}
