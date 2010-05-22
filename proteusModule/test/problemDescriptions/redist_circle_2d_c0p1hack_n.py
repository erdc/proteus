from pyadh import *
from pyadh.default_n import *
from redist_circle_2d_p import *

timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_A

runCFL = 0.1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=21
nLevels = 2
DT = None

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
