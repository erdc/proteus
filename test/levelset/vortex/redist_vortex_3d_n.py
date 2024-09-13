from proteus import *
from proteus.default_n import *
from .redist_vortex_3d_p import *
from .vortex import *

timeIntegration = NoIntegration
stepController = Newton_controller

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
        elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)    
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)
    if LevelModelType == RDLS.LevelModel:
        subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
    else:
        subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=False)
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=shockCapturingFactor_rd,lag=False)
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

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.0

nl_atol_res = atolRedistance

maxNonlinearIts = 100 #1 for PTC
maxLineSearches = 100
levelNonlinearSolverConvergenceTest='rits'

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    linear_solver_options_prefix = 'rdls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = {}
