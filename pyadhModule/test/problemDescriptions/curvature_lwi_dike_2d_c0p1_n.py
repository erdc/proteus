from pyadh import *
from pyadh.default_n import *
from curvature_lwi_dike_2d_p import *
from lwi_dike import *

timeIntegration = NoIntegration
timeIntegrator = SteadyStateIntegrator

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,lwi_dike_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,lwi_dike_quad_order)

subgridError = None
shockCapturing = None
massLumping = False
reactionLumping=False
numericalFluxType = Curvature_exterior
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False
tolFac = 0.0
nl_atol_res = 1.0e-8 #1e-4
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
linTolFac = 0.001
conservativeFlux = None
