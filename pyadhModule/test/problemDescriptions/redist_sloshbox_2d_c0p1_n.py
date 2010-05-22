from pyadh import *
from pyadh.default_n import *
from redist_sloshbox_2d_p import *
from sloshbox import *

if rdtimeIntegration == 'newton':    
    timeIntegration = NoIntegration
    stepController = Newton_controller
elif rdtimeIntegration == 'tte':
    timeIntegration = BackwardEuler_cfl
    timeIntegration = PsiTCtte
elif rdtimeIntegration == 'osher-fmm':
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_FMM_controller
    runCFL=1.0
else:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller
    runCFL=1.0#0.33
    rtol_res[0] = 0.0
    atol_res[0] = he*0.1#1.0e-6

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

if LevelModelType == RDLS.OneLevelRDLS and not RDLS.debugRDLS:
    subgridError.lag=0

#subgridError=None
if rdtimeIntegration == 'newton':    
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)
else:
    if LevelModelType == RDLS.OneLevelRDLS:
        shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)
    else:
        shockCapturing = Eikonal_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)
        shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)
        #shockCapturing = HamiltonJacobi_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)
        #shockCapturing = ConstantDiffusion_SC(coefficients,nd,shockCapturingFactor=0.01,lag=True)
    
#subgridError.lag=0

massLumping = False

if LevelModelType == RDLS.OneLevelRDLS:
    numericalFluxType = DoNothing
else:    
    numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind

multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver
multilevelNonlinearSolver  = Newton#LNI
levelNonlinearSolver = Newton
if rdtimeIntegration != 'newton':    
    maxLineSearches = 0
    levelNonlinearSolverConvergenceTest='rits'
else:
    maxLineSearches=100
    levelNonlinearSolverConvergenceTest='rits'
nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.0

if rdtimeIntegration != 'newton':
    nl_atol_res = rd_atol#0.05*he#1.0e-5
else:
    nl_atol_res = rd_atol# 0.05*he#1.0e-7#0.01*L[0]/nnx

maxNonlinearIts = 25#25 #1 for PTC

matrix = SparseMatrix

if LevelModelType == RDLS.OneLevelRDLS:
    numericalFluxType = DoNothing
else:    
    numericalFluxType = None

if usePETSc:
    numericalFluxType = DoNothing

    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
