from pyadh import *
from pyadh.default_n import *
from curvature_4petal_2d_p import *


timeIntegration = NoIntegration
timeIntegrator = ForwardIntegrator

runCFL = 2.0
DT = None
nDTout = 1


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
femSpaces = {0:DG_Constants}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)


nn=7
nLevels = 3

subgridError = None
shockCapturing = None
massLumping = False
numericalFluxType = None
numericalFluxType = Advection_Diagonal_average
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False
tolFac = 1.0
nl_atol_res = 1.0e-4 #1e-4
maxNonlinearIts = 100
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
linTolFac = 0.001
conservativeFlux = {0:'point-eval'}
