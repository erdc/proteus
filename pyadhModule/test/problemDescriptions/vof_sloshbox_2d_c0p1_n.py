from pyadh import *
from pyadh.default_n import *
from sloshbox import *
from vof_sloshbox_2d_p import *

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)


if useBackwardEuler_vof:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    if timeOrder == 2:
        timeIntegration = VBDF
        stepController = Min_dt_cfl_controller
        
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[0] = vof_rtol
    atol_u[0] = vof_atol
#timeIntegration = BackwardEuler_cfl
#stepController=Min_dt_controller
if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elif spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=vof_shockCapturingFactor,lag=True)#linear
subgridError = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)
massLumping = False
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 0.001*he#1.0e-8

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
   
linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
