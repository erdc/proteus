from pyadh import *
from pyadh.default_n import *
from dr_2c_di_het_ss_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:NC_AffineLinearOnSimplexWithNodalBasis}

#try to pick terms specifically for components and equation terms
# elementQuadrature = {'default':SimplexGaussQuadrature(nd,3),
#                      ('a',1,0):SimplexGaussQuadrature(nd,3), #pick midpoint quad for P1_NC components
#                      ('a',1,1):SimplexGaussQuadrature(nd,3),
#                      ('r',1):SimplexGaussQuadrature(nd,1),#try nd,2 to see how new pp works
#                      ('f',1):SimplexGaussQuadrature(nd,3)} #could change ('u',1):SimplexGaussQuadrature(nd,2) too

# elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3),
#                              ('a',1,0):SimplexGaussQuadrature(nd-1,3), #pick midpoint quad for P1_NC components
#                              ('a',1,1):SimplexGaussQuadrature(nd-1,3),
#                              ('r',1):SimplexGaussQuadrature(nd-1,1), #try nd-1,2 to see new pp effect
#                              ('f',1):SimplexGaussQuadrature(nd-1,3)} #could add ('u',1):SimplexGaussQuadrature(nd-1,1)
elementQuadrature = {'default':SimplexGaussQuadrature(nd,3)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3)}


nn = 5
nLevels = 4

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001
#{0:None,1:'p1-nc'}, {0:'pwl-bdm',1:'p1-nc'},{0:'pwl',1:'p1-nc'},{0:'point-eval',1:'p1-nc'},
conservativeFlux = {0:'pwl-bdm',1:'p1-nc'} 
