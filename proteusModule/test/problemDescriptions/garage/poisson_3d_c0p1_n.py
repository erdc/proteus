from pyadh import *
from pyadh.default_n import *
from poisson_3d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#try to pick terms specifically for components and equation terms
elementQuadrature = {'default':SimplexGaussQuadrature(nd,3)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3)}

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
elementQuadrature = SimplexGaussQuadrature(nd,2)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

# femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
# elementQuadrature = SimplexGaussQuadrature(nd,4)
# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

# femSpaces = {0:DG_AffineP4_OnSimplexWithMonomialBasis}
# elementQuadrature = SimplexGaussQuadrature(nd,5)
# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

#femSpaces = {0:DG_AffineP5_OnSimplexWithMonomialBasis}
#elementQuadrature = SimplexGaussQuadrature(nd,6)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,6)

numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#mwf everything gets SimplexGaussQuadrature
#elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)



nn = 3
nLevels = 2

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd)

shockCapturing = None

multilevelNonlinearSolver  = Newton#NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

atol = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
