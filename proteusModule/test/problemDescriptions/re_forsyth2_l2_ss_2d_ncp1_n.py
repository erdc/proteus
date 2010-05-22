from pyadh import *
from pyadh.default_n import *
from re_forsyth2_l2_ss_2d_p import *

timeIntegration = NoIntegration

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

#
elementQuadrature = {'default':SimplexGaussQuadrature(nd,2),
                     ('u',0):SimplexGaussQuadrature(nd,2)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,2),
                             ('u',0):SimplexGaussQuadrature(nd-1,2)}

nn=3
nLevels = 3

subgridError = None

massLumping = False

numericalFluxType = None

shockCapturing = None

multilevelNonlinearSolver = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

maxNonlinearIts = 1001

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'p1-nc'}
