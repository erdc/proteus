from pyadh import *
from pyadh.default_n import *
from gw_block_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

gw_quad_order = 3
#mwf everything gets SimplexGaussQuadrature
elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)

triangleOptions = "q30Dena0.005A"
nLevels = 1

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc

levelLinearSolver = LU
#levelLinearSolver = PETSc

linTolFac = 1.0e-10

#numericalFluxType = Diffusion_IIPG_exterior
conservativeFlux = {0:'pwl'}#{0:'pwc'}
