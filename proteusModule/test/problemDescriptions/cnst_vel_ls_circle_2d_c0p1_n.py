from pyadh import *
from pyadh.default_n import *
from cnst_vel_ls_circle_2d_p import *

timeIntegration = BackwardEuler

runCFL = 0.1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)


nn=21
nLevels = 2
DT = None
nDTout = 10


#subgridError = DiagonalStabilization_1
subgridError = HamiltonJacobi_ASGS(coefficients,nd)
#subgridError = None

#must also have basic quadrature set to SimplexLobattoQuadrature 
massLumping = False

numericalFluxType = None

shockCapturing = None

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
