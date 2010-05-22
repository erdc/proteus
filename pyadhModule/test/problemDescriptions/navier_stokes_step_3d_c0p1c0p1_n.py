from pyadh import *
from pyadh.default_n import *
from navier_stokes_step_3d_p import *

# timeIntegration = FLCBDF
# stepController  = FLCBDF_controller
# systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
timeIntegration = NoIntegration
stepController  = Newton_controller
timeIntegration = BackwardEuler
stepController  = FixedStep
# rtol_u[1] = 1.0e-2
# rtol_u[2] = 1.0e-2
# rtol_u[3] = 1.0e-2
# atol_u[1] = 1.0e-2
# atol_u[2] = 1.0e-2
# atol_u[3] = 1.0e-2

# timeIntegration = BackwardEuler_cfl
# stepController = Min_dt_controller
# systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
# runCFL = 1.0#
#DT=1.0e-1
#DT=0.0000001
nDTout = 1#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#             1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#             2:C0_AffineQuadraticOnSimplexWithNodalBasis,
#             3:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
triangleOptions = "VpAq1.25Dena%e" % ((0.5*upstream_height)**3 / 6.0)
print triangleOptions
#triangleOptions = "Aq1.5Den"
nLevels = 1

if LevelModelType == RANS2P.OneLevelRANS2P:
    subgridError = NavierStokesASGS_velocity_pressure_opt(coefficients,nd,lag=True,delayLagSteps=-1)
    shockCapturing = NavierStokes_SC_opt(coefficients,nd,0.25,lag=True)
else:
    subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=-1)
    shockCapturing = NavierStokes_SC(coefficients,nd,0.0,lag=True)

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior

massLumping = False


multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 1000

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc

levelLinearSolver = LU
#levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

if LevelModelType == RANS2P.OneLevelRANS2P:
    conservativeFlux = {0:'pwl-bdm-opt'}#None
else:
    conservativeFlux = {0:'pwl'}
#parallelPartitioningType = MeshParallelPartitioningTypes.element
#nLayersOfOverlapForParallel = 0
