from pyadh import *
from pyadh.default_n import *
from navier_stokes_manning_2d_tmp_p import *

timeIntegration = NoIntegration
stepController = Newton_controller
tnList=[0.0,1.0]
#timeIntegration = PsiTCtte_new
#stepController = PsiTCtte_controller
#systemStepControllerType = SplitOperator.Sequential_FixedStep
#timeIntegration = BackwardEuler_cfl
#stepController = Osher_controller
#tnList=[0.0,1.0e-4]
#runCFL = 0.33
#timeIntegration = FLCBDF
#stepController  = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_res[1] = 1.0e-2
rtol_res[2] = 1.0e-2
atol_res[1] = 1.0e-2
atol_res[2] = 1.0e-2
rtol_u[1] = 1.0e-2
rtol_u[2] = 1.0e-2
atol_u[1] = 1.0e-2
atol_u[2] = 1.0e-2
#tnList=[0.0,1.0e-2,1.0e6]





spaceOrder = 1
parallel = False
triangleOptions="pAq30Dena%f" % (0.0000025*(dx)**2,)

if spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 1:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 2:C0_AffineQuadraticOnSimplexWithNodalBasis}
    hFactor = 0.5
else:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis}
    hFactor = 1.0

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,3)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,3)

# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}
# elementQuadrature = SimplexGaussQuadrature(nd,5)

# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False,hFactor=hFactor)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd)

massLumping = False

#shockCapturing = NavierStokes_SC(coefficients,nd,0.5,lag=True)
#shockCapturing = None

matrix = SparseMatrix

if parallel:
    assert spaceOrder == 1, "P1 only in parallel for now"
    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    #default number of layers to use > 1 with element partition means
    #C0P1 methods don't need to do communication in global element assembly
    #nodal partitioning does not need communication for C0P1 (has overlap 1) regardless
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 2
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    #numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior
    numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU


multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
maxNonlinearIts =50
maxLineSearches =50
nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8#2.0e-4

auxiliaryVariables=[VelocityAverage()]
#auxiliaryVariables=[BoundaryPressure()]
