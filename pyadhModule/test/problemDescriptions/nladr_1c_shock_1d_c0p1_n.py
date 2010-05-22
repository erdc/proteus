from pyadh import *
from pyadh.default_n import *
from nladr_1c_shock_1d_p import *

timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_A
#timeIntegration = ForwardEuler

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 5
DT = None
nDTout = 100

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='1',lag=True)

massLumping = False

numericalFluxType = None

shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver =Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts =101

tolFac = 0.001

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = NI#LU

levelLinearSolver = MGM#LU

linearSmoother = GaussSeidel

linTolFac = 0.0001

conservativeFlux = None
