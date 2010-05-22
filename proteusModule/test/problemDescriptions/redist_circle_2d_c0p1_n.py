from pyadh import *
from pyadh.default_n import *
from redist_circle_2d_p import *

#timeIntegration = BackwardEuler
#timeIntegrator = ForwardIntegrator
#timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
#timeIntegrator  = SteadyStateIntegrator
#timeIntegrator  = ForwardIntegrator
#timeIntegration = BackwardEuler_cfl
timeIntegration = VBDF
timeOrder = 2
#stepController = Osher_controller
#mwf try combining Osher with FMM 
stepController = Osher_FMM_controller

nDTout = 1
DT = None
#small cfl don't freeze level set
runCFL = 0.01

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#elementQuadrature = SimplexLobattoQuadrature(nd,4)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,4)

nn=41#81
nLevels = 1


subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False,stabFlag='2')
if LevelModelType == RDLS.OneLevelRDLS:
        subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)

#works if freeze level set
#shockCapturing = ResGradQuad_SC(coefficients,nd,lag=False,shockCapturingFactor=0.99)
#works if don't freeze level set and take small steps
shockCapturing = ResGradQuad_SC(coefficients,nd,lag=False,shockCapturingFactor=0.2)

numericalFluxType = None
if LevelModelType == RDLS.OneLevelRDLS:
    numericalFluxType = DoNothing
multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#these need to be set appropriately for pseudo-transient
tolFac = 0.0
nl_atol_res = 1.0e-4
maxNonlinearIts =100

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
