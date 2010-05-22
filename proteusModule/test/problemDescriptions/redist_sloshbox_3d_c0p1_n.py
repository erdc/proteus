from pyadh import *
from pyadh.default_n import *
from redist_sloshbox_3d_p import *
from sloshbox3d import *

if rdtimeIntegration == 'newton':    
    timeIntegration = NoIntegration
    stepController = Newton_controller
else:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller
    runCFL=1.0
    rtol_res[0] = 0.0
    atol_res[0] = 0.1*he#1.0e-4

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
if spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

subgridErrorType = HamiltonJacobi_ASGS
if LevelModelType == RDLS.OneLevelRDLS and not RDLS.debugRDLS:
    subgridErrorType = HamiltonJacobi_ASGS_opt
if rdtimeIntegration == 'newton':
    subgridError = subgridErrorType(coefficients,nd,stabFlag='2',lag=False)
else:
    subgridError = subgridErrorType(coefficients,nd,stabFlag='2',lag=False)
subgridError.lag=0
if rdtimeIntegration == 'newton':    
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)
else:
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)
    
massLumping = False

if LevelModelType == RDLS.OneLevelRDLS:
    numericalFluxType = DoNothing
else:    
    numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind

#multilevelNonlinearSolver  = MultilevelEikonalSolver
#levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
if rdtimeIntegration != 'newton':    
    maxLineSearches = 0
    levelNonlinearSolverConvergenceTest='rits'
    maxNonlinearIts = 1 #1 for PTC
else:
    maxLineSearches=100
    levelNonlinearSolverConvergenceTest='rits'
    maxNonlinearIts = 25 #1 for PTC
    
nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.0

if rdtimeIntegration != 'newton':
    nl_atol_res = 0.01*he
else:
    nl_atol_res = 0.01*he#1.0e-4#0.01*L[0]/nnx

matrix = SparseMatrix

if usePETSc:
    numericalFluxType = DoNothing

    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    numericalFluxType = DoNothing

    multilevelLinearSolver = LU
    
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
