from pyadh import *
from pyadh.default_n import *
from ladr_ss_2d_p import *

timeIntegration = NoIntegration

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP5_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

# femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
# elementQuadrature = SimplexGaussQuadrature(nd,4)
# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
# elementQuadrature = SimplexGaussQuadrature(nd,4)
# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
#femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
#elementQuadrature = SimplexGaussQuadrature(nd,6)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,6)

# femSpaces = {0:DG_AffineP4_OnSimplexWithMonomialBasis}
# elementQuadrature = SimplexGaussQuadrature(nd,5)
# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

#femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
#elementQuadrature = SimplexGaussQuadrature(nd,6)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,6)

nn=5
nLevels = 1

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False,stabFlag='2')
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)
#shockCapturing = None
#massLumping = True


numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG_exterior
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG
#numericalFluxType = Advection_DiagonalUpwind

# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
# elementQuadrature = SimplexGaussQuadrature(nd,2)
# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)
# numericalFluxType = None


#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton
#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLJacobi
#multilevelNonlinearSolver  = NLGaussSeidel

#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver  = NLStarILU
#levelNonlinearSolver  = NLGaussSeidel
#levelNonlinearSolver  = NLJacobi

nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
#nonlinearSmoother = NLJacobi

fullNewtonFlag = True
#fullNewtonFlag = False

tolFac = 0.0#1.0e-8

nl_atol_res = 1.0e-8

maxNonlinearIts =1001

matrix = SparseMatrix
#matrix = Numeric.array

#multilevelLinearSolver = NI
multilevelLinearSolver = LU
#multilevelLinearSolver = StarILU
#multilevelLinearSolver = GaussSeidel
#multilevelLinearSolver = Jacobi
#multilevelLinearSolver = MGM

levelLinearSolver = LU
#levelLinearSolver = GaussSeidel
#levelLinearSolver = StarILU
#levelLinearSolver = MGM
#computeEigenvalues=True

#linearSmoother = StarILU#GaussSeidel
linearSmoother = GaussSeidel
#linearSmoother = Jacobi

linTolFac = 0.001

conservativeFlux = None
