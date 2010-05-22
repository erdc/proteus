from pyadh import *
from pyadh.default_n import *
from cnst_vel_ls_vortex_2d_p import *

#timeIntegration = BackwardEuler
timeIntegration = ForwardEuler_H
timeIntegration = OuterTheta
stepController = MinDT_controller

runCFL = 0.1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)


nn=3
nLevels = 5
DT = None
nDTout = 5

#subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd)

shockCapturing = None
#shockCaptuirng = ResGrad_SC(coefficients,nd)

massLumping = False

numericalFluxType = None

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
