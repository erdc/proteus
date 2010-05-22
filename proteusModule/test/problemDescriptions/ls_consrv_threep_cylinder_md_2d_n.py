from pyadh import *
from pyadh.default_n import *
from ls_consrv_threep_cylinder_md_2d_p import *


timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = None

massLumping = False

numericalFluxType = None

if usePETSc:
    numericalFluxType = Diffusion_IIPG_exterior

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 0.001*he#1.0e-10

maxNonlinearIts = 10

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    
linearSmoother = GaussSeidel

linTolFac = 1.0e-6

conservativeFlux = None
