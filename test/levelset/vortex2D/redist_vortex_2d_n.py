from proteus import *
from proteus.default_n import *
try:
    from .redist_vortex_2d_p import *
    from .vortex2D import *
except:
    from redist_vortex_2d_p import *
    from vortex2D import *

if redist_Newton:
    timeIntegration = NoIntegration
    stepController = Newton_controller
    tolFac = 0.0
    nl_atol_res = atolRedistance
    maxNonlinearIts = 100
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

if cDegree_ls==0:
    if useHex:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
        elementQuadrature = CubeGaussQuadrature(nd,vortex_quad_order)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,vortex_quad_order)
    else:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
        if compQuad:
            base_quad_rule = SimplexGaussQuadrature(nd,vortex_quad_order)
            elementQuadrature = CompositeTriangle(base_quad_rule,hk)
        else:
            elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)
    if LevelModelType == RDLS.LevelModel:
        subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
    else:
        subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=False)
    shockCapturing = RDLS.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_rd,lag=lag_shockCapturing_rd)
    if parallel or LevelModelType == RDLS.LevelModel:
        numericalFluxType = DoNothing
elif cDegree_ls==-1:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller
    runCFL=1.0
    rtol_res[0] = 0.0
    atol_res[0] = 0.01/float(nn-1)#1.0e-6
    if pDegree_ls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    numericalFluxType=HamiltonJacobi_DiagonalLesaintRaviart


multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None

fullNewtonFlag = True


matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py#PETSc
    levelLinearSolver = KSP_petsc4py#PETSc
    linear_solver_options_prefix = 'rdls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU


conservativeFlux = {}
