from proteus import *
from proteus.default_n import *
from moveMesh_p import *

timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)


nLevels = 1

subgridError = None

massLumping = False

numericalFluxType = Stress_IIPG_exterior

shockCapturing = None

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nonlinearSmoother = None
linearSmoother = None

fullNewtonFlag = True


matrix = SparseMatrix

if useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'mesh_'
linearSmoother = None
linearSolverConvergenceTest = 'r-true'
tolFac = 0.0
linTolFac = 0.001
l_atol_res = 0.001*ct.mesh_nl_atol_res
nl_atol_res = ct.mesh_nl_atol_res
maxNonlinearIts = 50#should be linear
maxLineSearches = 0

conservativeFlux = None

auxiliaryVariables=[fo]
