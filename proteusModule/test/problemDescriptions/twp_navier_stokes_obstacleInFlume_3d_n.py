from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_obstacleInFlume_3d_p import *
from obstacleInFlume3d import *

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
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
elementQuadrature = SimplexGaussQuadrature(nd,obstacleInFlume_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,obstacleInFlume_quad_order)

subgridError = None

#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=noPressureStabilization)
subgridError = NavierStokesASGS_velocity_pressure_opt(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=hFactor,noPressureStabilization=noPressureStabilization)

massLumping = False

#shockCapturing = NavierStokes_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)
shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=True)

numericalFluxType = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 25
maxLineSearches =50

tolFac = 0.0#0.01

nl_atol_res = 1.0e-4#0.001*he
linearSolverConvergenceTest = 'r-true'

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
