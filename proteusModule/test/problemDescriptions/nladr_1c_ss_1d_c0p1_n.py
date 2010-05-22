from pyadh import *
from pyadh.default_n import *
from nladr_1c_ss_1d_p import *

timeIntegration = NoIntegration

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nLevels = 7

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

numericalFluxType  = None
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLStarILU

fullNewtonFlag = True

tolFac = 0.01
maxNonlinearIts =1001

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
