from pyadh import *
from pyadh.default_n import *
from curvature_bubble_3d_p import *
from bubble3d import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
timeIntegrator = SteadyStateIntegrator
#DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,bubble_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
shockCapturing = None
massLumping = False
reactionLumping=False
numericalFluxType = None
#numericalFluxType = Advection_Diagonal_average
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False
tolFac = 0.0
atol = 1.0e-8 #1e-4
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
linTolFac = 0.001
conservativeFlux = None
