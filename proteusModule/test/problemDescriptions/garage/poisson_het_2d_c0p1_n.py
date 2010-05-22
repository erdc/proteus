from pyadh import *
from pyadh.default_n import *
from poisson_het_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#mwf everything gets SimplexGaussQuadrature
#cek makeing simple 1c problem
elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 5
nLevels = 2

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-8

atol = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

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

linTolFac = 1.0e-91

conservativeFlux = {0:'pwl',1:'p1-nc'} #{0:None,1:'p1-nc'}
