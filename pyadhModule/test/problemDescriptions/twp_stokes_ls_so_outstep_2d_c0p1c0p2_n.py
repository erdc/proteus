from pyadh import *
from pyadh.default_n import *
from twp_stokes_ls_so_outstep_2d_p import *

timeIntegration = NoIntegration

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineQuadraticOnSimplexWithNodalBasis,
             2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

nn=3
nLevels = 4

subgridError = None

massLumping = None

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 10

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
