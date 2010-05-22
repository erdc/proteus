from pyadh import *
from pyadh.default_n import *
from ls_vortex_2d_p import *
from vortex import *

timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-4
atol_u[0] = 1.0e-4

femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)


subgridError = None
if useHJ:
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)
else:
    subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=True)

massLumping = False

numericalFluxType = None

shockCapturing = None#ResGrad_SC(coefficients,nd)

multilevelNonlinearSolver  = Newton

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

#checkMass = True

