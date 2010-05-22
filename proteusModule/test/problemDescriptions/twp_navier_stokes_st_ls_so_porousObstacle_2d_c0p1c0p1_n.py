from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_st_ls_so_porousObstacle_2d_p import *
from porousObstacle import *

timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator

#DT=None
#runCFL=0.15
#DT=1.0e-5

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,porousObstacle_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,porousObstacle_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 10

tolFac = 0.1

nl_atol_res = 0.01/(nn-1.0)

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 1.0e-8

conservativeFlux = {0:'pwl',1:'point-eval',2:'point-eval'}
#conservativeFlux = {0:'point-eval',1:'point-eval',2:'point-eval'}
