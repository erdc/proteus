from pyadh import *
from pyadh.default_n import *
from redist_cnst_vel_circle_2d_p import *

timeIntegration = PsiTCtte
#timeIntegration = NoIntegration
timeIntegrator = SteadyStateIntegrator

DT = None
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=21
nLevels = 2

subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)

#shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)

massLumping = False

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.01
nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
