from pyadh import *
from pyadh.default_n import *
from curvature_wavetank_2d_p import *
from wavetank import *

timeIntegration = NoIntegration
stepController = FixedStep

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,wavetank_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,wavetank_quad_order)

#this is a linear problem
subgridError = None
shockCapturing = None
massLumping = False
reactionLumping=False
numericalFluxType = Curvature_exterior
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False
tolFac = 0.0
nl_atol_res = 1.0e-8
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = PETSc
levelLinearSolver = PETSc
linTolFac = 0.001
conservativeFlux = None
