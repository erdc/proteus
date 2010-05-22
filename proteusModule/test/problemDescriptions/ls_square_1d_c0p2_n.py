from pyadh import *
from pyadh.default_n import *
from ls_square_1d_p import *
from square import *

nn=26; nLevels = 1
#runCFL=0.3
timeIntegrator = ForwardIntegrator
#timeIntegration = BackwardEuler_cfl
#timeIntegration = ForwardEuler_A
#stepController  = Min_dt_controller
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-4
atol_u[0] = 1.0e-4

femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)


subgridError = Advection_ASGS(coefficients,nd,stabFlag='2')


massLumping = False

numericalFluxType = None

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 0.01/(nn -1.0 )

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {}



