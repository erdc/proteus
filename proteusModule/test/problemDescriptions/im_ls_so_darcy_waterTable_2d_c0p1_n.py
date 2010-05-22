from pyadh import *
from pyadh.default_n import *
from im_ls_so_darcy_waterTable_2d_p import *

timeIntegration = BackwardEuler
DT = None #0.1
nDTout = int(T/0.1)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 5

#subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = None

#shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

maxNonlinearIts = 100

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
