from pyadh import *
from pyadh.default_n import *
from la_shock_1d_p import *

timeIntegration = ForwardEuler_A
timeIntegration = ForwardEuler
stepController = Min_dt_controller
# timeIntegration = BackwardEuler
# stepController = FixedStep
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-2
atol_u[0] = 1.0e-2

runCFL = 0.01
DT = None
nDTout = 100
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 21
nLevels = 1

subgridError = Advection_ASGS(coefficients,nd,stabFlag='1',lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd)

multilevelNonlinearSolver  = Newton#NLNI

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

checkMass = True

