from pyadh import *
from pyadh.default_n import *
from gw_het_fivespot_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

#mwf everything gets SimplexGaussQuadrature
#cek makeing simple 1c problem
elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)
#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = PETSc
#levelLinearSolver = MGM
#levelLinearSolver = StarILU
#levelLinearSolver = GaussSeidel
#levelLinearSolver = Jacobi

smoother = StarILU
smoother = GaussSeidel
smoother = Jacobi

linTolFac = 1.0e-10

#numericalFluxType = Diffusion_IIPG_exterior
conservativeFlux = {0:'p1-nc'}#{0:'pwl'}
