from pyadh import *
from pyadh.default_n import *
from navier_stokes_squareobstacle_2d_p import *

timeIntegration = BackwardEuler
timeIntegration = FLCBDF
timeIntegrator = ForwardIntegrator
atol_u[0] = 1.0e-4
rtol_u[0] = 1.0e-4
timeOrder= 5

runCFL = 0.1#None
#DT=1.0e-1
DT=0.0000001
nDTout = 1#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 1.0e-5

nl_atol_res = 1.0e-5

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
