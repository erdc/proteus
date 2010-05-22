from pyadh import *
from pyadh.default_n import *
from twp_darcy_im_ls_so_waterTable_2d_p import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
#timeIntegrator = SteadyStateIntegrator
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#mwf everything gets SimplexGaussQuadrature
elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 5
nLevels = 2

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

#set for PsiTCtte
tolFac = 0.01
nl_atol_res = 1.0e-6
maxNonlinearIts =20

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = {0:'pwl-bdm'} 

