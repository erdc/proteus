from pyadh import *
from pyadh.default_n import *
from ls_wavetank_2d_p import *
from wavetank import *

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,wavetank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,wavetank_quad_order)

subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = NF_base

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.e-6

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
