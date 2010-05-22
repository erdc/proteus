from pyadh import *
from pyadh.default_n import *
from redist_simpleweir_2d_p import *
from simpleweir import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
timeIntegrator = SteadyStateIntegrator
DT=None
runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

subgridError = None
#subgridError = HamiltonJacobi_ASGS(coefficients,nd)
subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.0,lag=False)

massLumping = False

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 1.0
nl_atol_res = 1.0e-4
#maxNonlinearIts = 1 #1 for PTC
maxNonlinearIts = 100 #1 for PTC
matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
