from pyadh import *
from pyadh.default_n import *
from gw_het_utep_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#mwf everything gets SimplexGaussQuadrature
#cek makeing simple 1c problem
elementQuadrature = SimplexGaussQuadrature(nd,2)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)
#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

triangleOptions += "A"
nn = 3
nLevels = 1

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = MGM
#levelLinearSolver = StarILU
#levelLinearSolver = GaussSeidel
#levelLinearSolver = Jacobi

smoother = StarILU
smoother = GaussSeidel
smoother = Jacobi

linTolFac = 1.0e-10

numericalFluxType = Diffusion_IIPG_exterior
conservativeFlux = {0:'pwl'}
