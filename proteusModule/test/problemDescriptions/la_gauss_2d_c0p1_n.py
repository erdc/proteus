from pyadh import *
from pyadh.default_n import *
from la_gauss_2d_p import *

#timeIntegration = BackwardEuler
#timeIntegration = OuterTheta
timeIntegration = FLCBDF
timeIntegration = FLCBDF
stepController = FLCBDF_controller
atol_u[0] = 1.0e-3
rtol_u[0] = 1.0e-3
runCFL = 0.9
DT = None
nDTout = 20

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=41
nLevels = 1

subgridError = None
subgridError = Advection_ASGS(coefficients,nd)

numericalFluxType = None

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd)

multilevelNonlinearSolver  = NLNI

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
