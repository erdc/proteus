from pyadh import *
from pyadh.default_n import *
from transport_het_fivespot_2d_p import *

timeIntegration = BackwardEuler_cfl
DT = None
runCFL = 0.9
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#mwf everything gets SimplexGaussQuadrature
#cek makeing simple 1c problem
elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)
#elementQuadrature = SimplexLobattoQuadrature(nd,gw_quad_order)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,gw_quad_order)


nLevels = 1
massLumping=False
#subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

shockCapturing = ResGrad_SC(coefficients,nd,lag=True)

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton
fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU
#multilevelLinearSolver = NI
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

levelLinearSolver = LU
#levelLinearSolver = MGM
#levelLinearSolver = StarILU
#levelLinearSolver = GaussSeidel
#levelLinearSolver = Jacobi

smoother = StarILU
smoother = GaussSeidel
smoother = Jacobi

linTolFac = 1.0e-10

conservativeFlux = None
