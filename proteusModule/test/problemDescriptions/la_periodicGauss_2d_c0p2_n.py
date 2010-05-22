from pyadh import *
from pyadh.default_n import *
from la_periodicGauss_2d_p import *
from periodicGauss import *

timeIntegration = FLCBDF
stepController = FLCBDF_controller
atol_u[0] = 1.0e-4
rtol_u[0] = 1.0e-4
DT = None

femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,space_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,space_quad_order)


if useHJ:
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)
else:
    subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=True)

numericalFluxType = None

if useHJ:
    shockCapturing = None#HamiltonJacobi_SC(coefficients,nd,lag=True)#None
else:
    shockCapturing = None#ResGrad_SC(coefficients,nd,lag=True)#None
    

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8
maxNonlinearIts =1

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True
