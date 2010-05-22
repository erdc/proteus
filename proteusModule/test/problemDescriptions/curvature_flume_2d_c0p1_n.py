from pyadh import *
from pyadh.default_n import *
from curvature_flume_2d_p import *
from flume import *

timeIntegration = NoIntegration
timeIntegrator = SteadyStateIntegrator

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,flume_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,flume_quad_order)

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
multilevelLinearSolver = PETSc
levelLinearSolver = PETSc
linTolFac = 0.001
conservativeFlux = None
