from pyadh import *
from pyadh.default_n import *
from ladr_2c_lw_ss_1d_p import *

timeIntegration = NoIntegration

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nLevels = 5

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
