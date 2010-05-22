from pyadh import *
from pyadh.default_n import *
from redist_lwi_dike_2d_p import *
from lwi_dike import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
timeIntegrator = SteadyStateIntegrator
timeIntegration = BackwardEuler
timeIntegrator = SignedDistanceIntegrator
#runCFL=1.0
#DT=None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,lwi_dike_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,lwi_dike_quad_order)

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

nl_atol_res = 0.01

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
