from pyadh import *
from pyadh.default_n import *
from curvature_threep_cylinder_md_2d_p import *
from threep_cylinder_md_2d import *

timeIntegration = NoIntegration
timeIntegrator = SteadyStateIntegrator

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

subgridError = None
shockCapturing = None
massLumping = False
reactionLumping=False
numericalFluxType = None
numericalFluxType = Curvature_exterior
multilevelNonlinearSolver  = Newton
fullNewtonFlag = True
tolFac = 0.0
nl_atol_res = 1.0e-8
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc

linTolFac = 0.001
conservativeFlux = None
