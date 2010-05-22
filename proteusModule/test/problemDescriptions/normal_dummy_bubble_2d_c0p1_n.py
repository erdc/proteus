from pyadh import *
from pyadh.default_n import *
from normal_dummy_bubble_2d_p import *
from bubble import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
timeIntegrator = SteadyStateIntegrator
DT=T/10.0
nDTout=10
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)


subgridError = None
shockCapturing = None
massLumping = False
#reactionLumping=True
numericalFluxType = Diffusion_IIPG_exterior
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False
tolFac = 1.0e-10
nl_atol_res = 1.0e-10 #1e-4
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
linTolFac = 0.001
conservativeFlux = None
