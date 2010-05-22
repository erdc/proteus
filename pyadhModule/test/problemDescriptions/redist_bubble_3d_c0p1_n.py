from pyadh import *
from pyadh.default_n import *
from redist_bubble_3d_p import *
from bubble3d import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
timeIntegrator = SteadyStateIntegrator
timeIntegration = BackwardEuler
timeIntegrator = SignedDistanceIntegrator
#runCFL=1.0
#DT=None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,bubble_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=False)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)

massLumping = False

numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind

#multilevelNonlinearSolver  = MultilevelEikonalSolver
#levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver
multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.01

nl_atol_res = 0.1*L[0]**3

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
