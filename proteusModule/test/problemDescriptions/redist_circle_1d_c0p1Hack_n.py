from pyadh import *
from pyadh.default_n import *
from redist_circle_1d_p import *

#timeIntegration = BackwardEuler
timeIntegration = ForwardEuler_A

runCFL = 0.1
DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=3
nLevels = 5

subgridError = HamiltonJacobi_ASGS(coefficients,nd)

massLumping = False

shockCapturing = ResGrad_SC(coefficients,nd)

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
