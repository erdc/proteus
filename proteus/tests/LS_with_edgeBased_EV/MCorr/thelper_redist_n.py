from proteus import *
from proteus.default_n import *
from thelper_redist_p import *
from thelper_cons_ls import *

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
fullNewtonFlag = True

if redist_Newton:
    timeIntegration = NoIntegration
    stepController = Newton_controller
    tolFac = 0.0
    nl_atol_res = atolRedistance
    maxNonlinearIts = 25
    maxLineSearches = 0
    useEisenstatWalker = True
    linTolFac = 0.0
    levelNonlinearSolverConvergenceTest='r'
    nonlinearSolverConvergenceTest='r'
else:
    timeIntegration = BackwardEuler_cfl
    stepController = RDLS.PsiTC
    runCFL=0.33
    psitc['nStepsForce']=3
    psitc['nStepsMax']=15
    psitc['reduceRatio']=3.0
    psitc['startRatio']=1.
    tolFac = 0.0
    rtol_res[0] = 0.0
    atol_res[0] = atolRedistance
    nl_atol_res = 0.01*atolRedistance
    maxNonlinearIts = 1
    maxLineSearches = 0
    useEisenstatWalker = True
    linTolFac = 0.0
    levelNonlinearSolverConvergenceTest='rits'
    nonlinearSolverConvergenceTest='rits'

if useHex:
    if pDegree_ncls==1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_ncls==2:
        femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    elementQuadrature = CubeGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,quad_order)
else:    
    if pDegree_ncls==1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ncls==2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    elementQuadrature = SimplexGaussQuadrature(nd,quad_order)    
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
shockCapturing = RDLS.ShockCapturing(coefficients,
                                     nd,
                                     shockCapturingFactor=shockCapturingFactor_rd,
                                     lag=lag_shockCapturing_rd)
if parallel or LevelModelType == RDLS.LevelModel:
    numericalFluxType = DoNothing

nonlinearSmoother = None
matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'rdls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

conservativeFlux = {}
