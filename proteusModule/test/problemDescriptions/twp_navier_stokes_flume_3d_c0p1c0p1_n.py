from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_flume_3d_p import *
from flume3d import *

timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,flume_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,flume_quad_order)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 20

tolFac = 0.01

nl_atol_res = 0.01

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 1.0e-8

conservativeFlux = {0:'pwl'}
