from pyadh import *
from pyadh.default_n import *
from stokes_duct_ss_3d_p import *

timeIntegration = NoIntegration
timeIntegrator = SteadyStateIntegrator

runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 11
nLevels = 1

subgridError = StokesASGS_velocity(coefficients,nd)
#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd)

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-5

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None
