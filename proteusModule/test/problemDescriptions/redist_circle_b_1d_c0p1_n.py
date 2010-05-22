from pyadh import *
from pyadh.default_n import *
from redist_circle_b_1d_p import *

#timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_H
timeIntegration = PsiTCtte
timeIntegrator = SteadyStateIntegrator

runCFL = 0.1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 5
DT = None

subgridError = HamiltonJacobi_ASGS(coefficients,nd)

massLumping = False

numericalFluxType = None

shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#these need to be set for pseudo transient tol
tolFac = 0.1
nl_atol_res = 1.0e-4

maxNonlinearIts = 1000

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
