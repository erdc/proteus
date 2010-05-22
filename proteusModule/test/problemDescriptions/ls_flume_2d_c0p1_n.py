from pyadh import *
from pyadh.default_n import *
from flume import *
from ls_flume_2d_p import *

timeIntegration = BackwardEuler_cfl
  
stepController = Min_dt_controller
#timeIntegration = FLCBDF
#stepController  = FLCBDF_controller
#rtol_u[0] = 1.0e-2
#atol_u[0] = 1.0e-2

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,flume_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,flume_quad_order)

subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=ls_shockCapturingFactor,lag=lag_ls_shockCapturing)

numericalFluxType = NF_base # does nothing

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
