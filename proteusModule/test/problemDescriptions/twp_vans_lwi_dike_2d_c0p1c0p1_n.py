from pyadh import *
from pyadh.default_n import *
from twp_vans_lwi_dike_2d_p import *
from lwi_dike import *

timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,lwi_dike_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,lwi_dike_quad_order)

if useStokes:
    subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)
else:
    subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 10

tolFac = 0.01

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 1.0e-8

conservativeFlux = {0:'pwl'}
