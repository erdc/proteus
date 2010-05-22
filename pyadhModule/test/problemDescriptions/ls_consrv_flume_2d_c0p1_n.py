from pyadh import *
from pyadh.default_n import *
from flume import *
from ls_consrv_flume_2d_p import *


timeIntegration = BackwardEuler#NoIntegration
timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,flume_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,flume_quad_order)

# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = Diffusion_IIPG_exterior

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0
nl_atol_res = 1.0e-8

maxNonlinearIts = 100

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 1.0e-6

conservativeFlux = None
