from pyadh import *
from pyadh.default_n import *
from ls_so_simpleweir_2d_p import *

timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_H
timeIntegrator = ForwardIntegrator
#DT=None
runCFL=None#2.0

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)


subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = None
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.1

nl_atol_res = 1.0e-4

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
