from pyadh import *
from pyadh.default_n import *
from poisson_2c_het_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:NC_AffineLinearOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
#tell coefficients class to use l2 projection on p1nc component
l2proj = {0:False,1:True}
coefficients.l2proj  = l2proj

nn = 5
nLevels = 5


subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

atol = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001
#{0:None,1:'p1-nc'}, {0:'pwl-bdm',1:'p1-nc'},{0:'pwl',1:'p1-nc'},{0:'point-eval',1:'p1-nc'},
conservativeFlux = {0:'pwl-bdm',1:'p1-nc'} 
