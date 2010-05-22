from pyadh import *
from pyadh.default_n import *
from twp_stokes_ls_fc_weir_2d_p import *

#timeIntegration = BackwardEuler
timeIntegration = ForwardEuler_A

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineQuadraticOnSimplexWithNodalBasis,
             3:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementQuadrature = SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,3)

nn = 21
polyfile = "weir"
nLevels = 1
DT = None

subgridError =TwophaseStokes_LS_FC_ASGS(coefficients,nd)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
