from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_cylinder_md_2d_p import *

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#             1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#             2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

if useStokes:
 subgridError = StokesASGS_velocity_pressure(coefficients,nd)
else:
 subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=0,hFactor=1.0)
#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=0,hFactor=1.0)

#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

massLumping = False

shockCapturing = NavierStokes_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 25

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU
levelLinearSolver = LU

#multilevelLinearSolver = PETSc
#levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'pwl'}#,1:'pwl-bdm',2:'pwl-bdm'}
#conservativeFlux = {0:'pwc'}#,1:'pwc',2:'pwc'}
#conservativeFlux = {0:'point-eval',1:'point-eval',2:'point-eval'}
auxiliaryVariables=[rc]
