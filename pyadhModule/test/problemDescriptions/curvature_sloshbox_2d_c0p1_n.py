from pyadh import *
from pyadh.default_n import *
from curvature_sloshbox_2d_p import *
from sloshbox import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
timeIntegrator = SteadyStateIntegrator
#DT = None

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
else:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    
#femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
shockCapturing = None
massLumping = False
reactionLumping=False
numericalFluxType = None
#numericalFluxType = Advection_Diagonal_average
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
numericalFluxType = Curvature_exterior
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False
tolFac = 0.0
nl_atol_res = 0.001*he#1.0e-8 #1e-4
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc

linTolFac = 0.001
conservativeFlux = None
