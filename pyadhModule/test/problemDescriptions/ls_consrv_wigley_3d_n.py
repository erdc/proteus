from pyadh import *
from pyadh.default_n import *
from ls_consrv_wigley_3d_p import *


timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = None

massLumping = False

numericalFluxType = DoNothing

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 0.001*he#1.0e-10

maxNonlinearIts = 10

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linearSmoother = StarILU
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 1.0e-6

conservativeFlux = None
