from pyadh import *
from pyadh.default_n import *
from wigley import *
from redist_wigley_3d_p import *

if rdtimeIntegration == 'newton':    
    timeIntegration = NoIntegration
    stepController = Newton_controller
else:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller
    runCFL=1.0
    rtol_res[0] = 0.0
    atol_res[0] = he*0.1#1.0e-6

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
if spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

if rdtimeIntegration != 'newton':    
    subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=True)
    
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)
else:
    subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)

massLumping = False

numericalFluxType = DoNothing

multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver
multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton
if rdtimeIntegration != 'newton':    
    maxLineSearches = 0
nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.0

if rdtimeIntegration != 'newton':
    nl_atol_res = 0.1*he
else:
    nl_atol_res = 0.1*he#1.0e-7#0.01*L[0]/nnx

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linearSmoother = StarILU
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
