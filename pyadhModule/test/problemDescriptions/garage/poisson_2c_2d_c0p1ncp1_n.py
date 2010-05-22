from pyadh import *
from pyadh.default_n import *
from poisson_2c_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:NC_AffineLinearOnSimplexWithNodalBasis}
femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis,1:NC_AffineLinearOnSimplexWithNodalBasis}
numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#mwf everything gets SimplexGaussQuadrature
#elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#try to pick terms specifically for components and equation terms
elementQuadrature = {'default':SimplexGaussQuadrature(nd,3),
                     ('a',1,0):SimplexGaussQuadrature(nd,1), #pick midpoint quad for P1_NC components
                     ('a',1,1):SimplexGaussQuadrature(nd,1),
                     ('r',1):SimplexGaussQuadrature(nd,1),
                     ('f',1):SimplexGaussQuadrature(nd,1)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3),
                             ('a',1,0):SimplexGaussQuadrature(nd-1,1), #pick midpoint quad for P1_NC components
                             ('a',1,1):SimplexGaussQuadrature(nd-1,1),
                             ('r',1):SimplexGaussQuadrature(nd-1,1),
                             ('f',1):SimplexGaussQuadrature(nd-1,1)}


nn = 11
nLevels = 2

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton#NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

atol = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = {0:'pwl',1:'p1-nc'}
conservativeFlux = {0:'pwl',1:'p1-nc'}
