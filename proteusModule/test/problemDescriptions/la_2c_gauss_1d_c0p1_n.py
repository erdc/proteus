from pyadh import *
from pyadh.default_n import *
from la_2c_gauss_1d_p import *

#BackwardEuler
timeIntegration = BackwardEuler
#timeIntegration = FLCBDF
#timeIntegration = ForwardEuler_A
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-4
rtol_u[1] = 1.0e-4
atol_u[0] = 1.0e-4
atol_u[1] = 1.0e-4


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

DT=None
runCFL = 0.1
#mwf tolerances for FLCBDF
atol_u[0]=1.0e-4#1.0e-6
rtol_u[0]=1.0e-4#1.0e-6

nDTout = 1

nn = 101
nLevels = 1

subgridError = None
subgridError = Advection_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.01,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-2

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {}



