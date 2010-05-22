from pyadh import *
from pyadh.default_n import *
from wavetank3d import *
from vof_wavetank_3d_p import *


timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,wavetank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,wavetank_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
subgridError = Advection_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

shockCapturing = None

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver= PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
