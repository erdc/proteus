from pyadh import *
from pyadh.default_n import *
from navier_stokes_tube_ss_2d_p import *

timeIntegration = NoIntegration
timeIntegrator=SteadyStateIntegrator

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nny=41
nnx=161
nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

maxNonlinearIts=25

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
