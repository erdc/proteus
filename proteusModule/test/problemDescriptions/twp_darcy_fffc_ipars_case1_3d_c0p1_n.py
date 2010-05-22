from pyadh import *
from pyadh.default_n import *
from twp_darcy_fffc_ipars_case1_3d_p import *

runCFL=None
timeIntegration = FLCBDF
stepController = FLCBDF_controller

rtol_u[0] = 1.0e-4
atol_u[0] = 1.0e-4
rtol_u[1] = 1.0e-4
atol_u[1] = 1.0e-4

tnList=[0.0,T]

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

elementQuadrature = SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

triangleOptions = "pAYfen"

nLevels = 1

massLumping=False

#subgridError = FFDarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=False)

subgridError = None

shockCapturing=None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 20#025

maxLineSearches = 10#0

fullNewtonFlag = True

tolFac = 1.0e-3

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

#conservativeFlux = {0:'pwl',1:'pwl'}
