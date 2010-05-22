from pyadh import *
from pyadh.default_n import *
from twp_navier_stokes_ls_so_susan_2d_p import *

timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

maxNonlinearIts = 5

tolFac = 0.01

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
