from proteus import *
from proteus.default_n import *
try:
    from .ncls_p import *
    from .vortex2D import *
except:
    from ncls_p import *
    from vortex2D import *
nd = 2

# About time integration 
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

# About nonlinear solver
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nl_atol_res = atolLevelSet
l_atol_res = atolLevelSet
tolFac = 0.0

maxNonlinearIts = 25
maxLineSearches = 0
fullNewtonFlag = True

if useHex:
    hex=True
    if pDegree_ls == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_ls == 2:
        femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print("pDegree_ls = %s not recognized " % pDegree_ls)
    elementQuadrature = CubeGaussQuadrature(nd,vortex_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,vortex_quad_order)
else:
    if pDegree_ls == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ls == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print("pDegree_ls = %s not recognized " % pDegree_ls)
    elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=True)
massLumping = False
numericalFluxType = None

shockCapturing = NCLS.ShockCapturing(coefficients,
                                     nd,
                                     shockCapturingFactor=shockCapturingFactor_ls,
                                     lag=lag_shockCapturing_ls)
numericalFluxType = DoNothing

nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
nonlinearSmoother = None

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py#PETSc
    levelLinearSolver = KSP_petsc4py#PETSc
    linear_solver_options_prefix = 'ncls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

conservativeFlux = {}
