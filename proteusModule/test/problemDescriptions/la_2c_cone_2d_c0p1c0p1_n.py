from pyadh import *
from pyadh.default_n import *
from la_2c_cone_2d_p import *

#timeIntegration = BackwardEuler
timeIntegration = ForwardEuler_A

runCFL = 0.1
DT = None
nDTout = 10
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 4

subgridError = None
subgridError = Advection_ASGS(coefficients,nd)

massLumping = False

numericalFluxType = None

shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd)

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

checkMass = False
