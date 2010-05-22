from pyadh import *
from pyadh.default_n import *
from cnst_vel_ls_rotatingcircle_2d_p import *

timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_A

runCFL = 0.1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=3
nLevels = 5

DT = None
nDTout = 20


subgridError = HamiltonJacobi_ASGS(coefficients,nd)

shockCapturing = ResGrad_SC(coefficients,nd)

massLumping = False

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
