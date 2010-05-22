from pyadh import *
from pyadh.default_n import *
from threep_cylinder_md_2d import *
from vof_threep_cylinder_md_2d_p import *

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)


timeIntegration = BackwardEuler
#stepController=Min_dt_controller
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

multilevelLinearSolver = LU

levelLinearSolver = LU

if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
   
linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
