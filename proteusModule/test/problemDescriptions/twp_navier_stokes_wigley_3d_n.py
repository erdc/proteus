from pyadh import *
from pyadh.default_n import *
from wigley import *
from twp_navier_stokes_wigley_3d_p import *

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
stepController = HeuristicNL_dt_controller
nonlinearIterationsFloor = 2
nonlinearIterationsCeil=4
nonlinearIterationsFloor = 2
nonlinearIterationsCeil=4
dtNLgrowFactor  = 1.5
dtNLreduceFactor= 0.5#75
                     
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = NavierStokesASGS_velocity_pressure_opt(coefficients,nd,lag=lag_ns_subgridError,delayLagSteps=1,hFactor=1.0)
#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=1,hFactor=1.0)

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior #need weak for parallel and global conservation

massLumping = False

shockCapturing = NavierStokes_SC_opt(coefficients,nd,ns_shockCapturingFactor,lag=lag_ns_shockCapturing)
#shockCapturing = NavierStokes_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 10
maxLineSearches =0

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-4#0.0001*he

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linearSmoother=StarILU
    linearSmoother=SimpleNavierStokes3D
    linear_solver_options_prefix = 'navierStokes_'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU


linTolFac = 0.001

conservativeFlux = {0:'pwl-bdm-opt'}
#conservativeFlux = {0:'pwl'}#,1:'pwl-bdm',2:'pwl-bdm'}
#conservativeFlux = {0:'pwc'}#,1:'pwc',2:'pwc'}
#conservativeFlux = {0:'point-eval',1:'point-eval',2:'point-eval'}
#auxiliaryVariables=[rc]
