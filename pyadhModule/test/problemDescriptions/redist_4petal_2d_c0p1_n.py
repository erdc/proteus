from pyadh import *
from pyadh.default_n import *
from redist_4petal_2d_p import *




#timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
#timeIntegrator = SteadyStateIntegrator
timeIntegration = BackwardEuler_cfl
#stepController = Osher_controller
#mwf try combining Osher with FMM 
stepController = Osher_FMM_controller

runCFL = 0.001#2.0
DT = None
nDTout = 1


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=51#7
nLevels = 1#5

subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False,stabFlag='2')
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='1')


shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=1.5,lag=True)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.2,lag=False)
#shockCapturing = HamiltonJacobi_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)

massLumping = False

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 1.0

nl_atol_res = 1.0e-4 #1e-4
maxNonlinearIts = 100

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
