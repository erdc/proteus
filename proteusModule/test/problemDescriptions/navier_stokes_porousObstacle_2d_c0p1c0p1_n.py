from pyadh import *
from pyadh.default_n import *
from navier_stokes_porousObstacle_2d_p import *

timeIntegration = BackwardEuler

runCFL = 10.0

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn=41
nLevels = 1
DT=1.0e-3#1.0e-3#1.0e-1#
nDTout = int(T/DT)#int(T/DT)



subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd)


shockCapturing = None

maxNonlinearIts =100

multilevelNonlinearSolver  = Newton#NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = {0:'pwl',1:'point-eval',2:'point-eval'}

