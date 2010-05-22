from pyadh import *
from pyadh.default_n import *
from navier_stokes_sphere_3d_p import *

timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
#timeIntegration = BackwardEuler
#stepController  = FixedStep
#systemStepControllerType = SplitOperator.Sequential_FixedStep
#timeIntegration = NoIntegration
#stepController  = Newton_controller
rtol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
rtol_u[3] = 1.0e-3
atol_u[1] = 1.0e-6
atol_u[2] = 1.0e-6
atol_u[3] = 1.0e-6
#runCFL=0.5

#DT=1.0e-1
#nDTout=int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 1
triangleOptions="VpAq1.25ena%e" % (((0.1*L[2])**3)/6.0,)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True) #None

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts = 50

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

#conservativeFlux = {0:'pwl',1:'point-eval',2:'point-eval',3:'point-eval'}
