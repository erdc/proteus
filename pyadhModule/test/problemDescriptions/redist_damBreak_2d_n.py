from pyadh import *
from pyadh.default_n import *
from redist_damBreak_2d_p import *
from damBreak import *

if rdtimeIntegration == 'newton':    
    timeIntegration = NoIntegration
    stepController = Newton_controller
else:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller
    runCFL=1.0
    rtol_res[0] = 0.0
    atol_res[0] = he*0.01#1.0e-6

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
if spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,damBreak_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,damBreak_quad_order)

if rdtimeIntegration == 'newton':
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=False)
else:
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=True)

if rdtimeIntegration == 'newton':    
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)
else:
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)
    
massLumping = False

numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind

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
    nl_atol_res = 1.0e-5
else:
    nl_atol_res = 0.001*he#1.0e-7#0.01*L[0]/nnx

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

if usePETSc:
    numericalFluxType = NF_base# does nothing

    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
