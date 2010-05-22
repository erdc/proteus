from pyadh import *
from pyadh.default_n import *
from poisson_2c_het_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:C0_AffineLinearOnSimplexWithNodalBasis}


elementQuadrature = {'default':SimplexGaussQuadrature(nd,4)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,4)}

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
conservativeFlux = {0:'pwc',1:'sun-rt0'}
#{0:'pwc',1:'pwl'}
#{0:'point-eval',1:'pwl'} 
#{0:'pwl-bdm',1:'pwl'} 
