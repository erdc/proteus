from pyadh import *
from pyadh.default_n import *
from curvature_bubble_2d_p import *
from bubble import *

timeIntegration = NoIntegration
stepController = Newton_controller
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,bubble_quad_order)
subgridError = None
shockCapturing = None
massLumping = False
reactionLumping=False
numericalFluxType = None
numericalFluxType = Curvature_exterior
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False
tolFac = 0.0
atol = 1.0e-8
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
linTolFac = 0.001
conservativeFlux = None
needEBQ_GLOBAL = True
if usePETSc:
    numericalFluxType = Curvature_exterior

    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
