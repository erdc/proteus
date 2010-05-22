from pyadh import *
from pyadh.default_n import *
from bubble import *
from ls_consrv_bubble_2d_p import *


timeIntegration = NoIntegration
stepController = Newton_controller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,bubble_quad_order)

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

nl_atol_res = 1.0e-10

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
