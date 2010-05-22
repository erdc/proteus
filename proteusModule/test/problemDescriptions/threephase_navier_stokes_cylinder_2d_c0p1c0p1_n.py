from pyadh import *
from pyadh.default_n import *
from threephase_navier_stokes_cylinder_2d_p import *

timeIntegration = NoIntegration
timeIntegrator=SteadyStateIntegrator
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
#timeIntegration = BackwardEuler
#stepController  = FixedStep
#systemStepControllerType = SplitOperator.Sequential_FixedStep
# timeIntegration = NoIntegration
# stepController  = Newton_controller
rtol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
atol_u[1] = 1.0e-2
atol_u[2] = 1.0e-2

nDTout=100
DT = T/float(nDTout)
tnList = [i*DT for i in range(nDTout+1)]

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nny=41
nnx=161
nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

maxNonlinearIts=5

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
