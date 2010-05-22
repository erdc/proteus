from pyadh import *
from pyadh.default_n import *
from poisson_2c_het_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:NC_AffineLinearOnSimplexWithNodalBasis}


#try to pick terms specifically for components and equation terms
elementQuadrature = {'default':SimplexGaussQuadrature(nd,3),
                     ('a',1,0):SimplexGaussQuadrature(nd,2), #pick midpoint quad for P1_NC components
                     ('a',1,1):SimplexGaussQuadrature(nd,2),
                     ('r',1):SimplexGaussQuadrature(nd,1),#can raise order for p1-nc v2 post processing
                     ('f',1):SimplexGaussQuadrature(nd,2)} #could change ('u',1):SimplexGaussQuadrature(nd,2) too

elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3),
                             ('a',1,0):SimplexGaussQuadrature(nd-1,2), #pick midpoint quad for P1_NC components
                             ('a',1,1):SimplexGaussQuadrature(nd-1,2),
                             ('r',1):SimplexGaussQuadrature(nd-1,1),#can raise order for p1-nc v2 post processing
                             ('f',1):SimplexGaussQuadrature(nd-1,2)} #could add ('u',1):SimplexGaussQuadrature(nd-1,1)

#elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 5
nLevels = 3


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
