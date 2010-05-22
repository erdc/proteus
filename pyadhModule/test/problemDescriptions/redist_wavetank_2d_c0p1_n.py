from pyadh import *
from pyadh.default_n import *
from redist_wavetank_2d_p import *
from wavetank import *

timeIntegration = BackwardEuler_cfl
stepController = Osher_PsiTC_controller
runCFL=1.0
rtol_res[0] = 0.0
atol_res[0] = 0.01*he#1.0e-4

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,wavetank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,wavetank_quad_order)

subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True,stabFlag='2')
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False,stabFlag='2')

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)

massLumping = False

numericalFluxType = NF_base

multilevelNonlinearSolver  = NLNI #Newton for PTC

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.0

nl_atol_res = 1.0e-5

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
