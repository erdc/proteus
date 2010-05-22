from pyadh import *
from pyadh.default_n import *
from linear_elasticity_plate_with_hole_2d_p import *

timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nLevels = 1
he=0.05
triangleOptions="pAq30Dena%f" % (0.5*he**2,)

subgridError = None

massLumping = False

numericalFluxType = None

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
