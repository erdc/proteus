from pyadh import *
from pyadh.default_n import *
from poisson_2c_ss_1d_p import *

timeIntegration = NoIntegration

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:NC_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 5

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

atol = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
