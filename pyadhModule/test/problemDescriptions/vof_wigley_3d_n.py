from pyadh import *
from pyadh.default_n import *
from wigley import *
from vof_wigley_3d_p import *

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)


timeIntegration = BackwardEuler_cfl
stepController=Min_dt_controller
stepController = HeuristicNL_dt_controller
nonlinearIterationsFloor = 3
nonlinearIterationsCeil=5
dtNLgrowFactor  = 1.5
dtNLreduceFactor= 0.5#75
#stepController = HeuristicNL_dt_controller
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
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

if usePETSc:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linearSmoother = StarILU
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
