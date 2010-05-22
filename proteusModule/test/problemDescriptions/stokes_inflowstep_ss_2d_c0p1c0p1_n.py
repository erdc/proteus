from pyadh import *
from pyadh.default_n import *
from stokes_inflowstep_ss_2d_p import *

timeIntegration = NoIntegration

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

nn=3
nLevels = 4

subgridError = StokesASGS_velocity(coefficients,nd)

massLumping = None

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 10

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = Numeric.array#SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
