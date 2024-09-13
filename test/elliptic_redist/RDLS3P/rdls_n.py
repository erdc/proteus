from proteus import *
from proteus.default_n import *
from .rdls_p import *
from .vortex2D import *

timeIntegration = NoIntegration
stepController = Newton_controller

# About the nonlinear solver
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = TwoStageNewton

tolFac = 0.0
nl_atol_res = atolRedistance
linTolFac = 0.0

maxNonlinearIts = 100000
maxLineSearches = 0
useEisenstatWalker = True
fullNewtonFlag = True

if useHex:
    hex=True
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

subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=False)
shockCapturing = RDLS.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_rd,lag=lag_shockCapturing_rd)

numericalFluxType = DoNothing

nonlinearSmoother = None
levelNonlinearSolverConvergenceTest='r'
nonlinearSolverConvergenceTest='r'

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
