from pyadh import *
from pyadh.default_n import *
from redist_flume_2d_p import *
from flume import *

timeIntegration = BackwardEuler_cfl
stepController = Osher_PsiTC_controller
runCFL=0.5
rtol_res[0] = 0.0
atol_res[0] = 1.0e-5

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,flume_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,flume_quad_order)

subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=True)

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)

numericalFluxType = NF_base

multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.0

nl_atol_res = 1.0e-6

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
