from pyadh import *
from pyadh.default_n import *
from sw_beach_erosion_board_1d_p import *

parallel = False
spaceOrder = 1; lagrangeBasis=True#False
if spaceOrder == 1:
    runCFL=0.1#0.3
elif spaceOrder == 2:
    runCFL=0.1
elif spaceOrder == 3:
    runCFL=0.05
else:
    runCFL=0.15
#runCFL=0.15
timeOrder = min(spaceOrder+1,3)
nStagesTime=timeOrder
timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller



if not lagrangeBasis and spaceOrder > 0:
    #limiterType = DGlimiterPkMonomial1d
    limiterType = DGlimiterPkMonomial1d_Sw
elif spaceOrder == 1:
    limiterType = DGlimiterP1Lagrange1d
    #limiterType = DGlimiterP1Lagrange1d_Sw
elif spaceOrder == 2:
    #limiterType = DGlimiterP2Lagrange1d
    limiterType = DGlimiterP2Lagrange1d_Sw

systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

if spaceOrder == 1:
    if lagrangeBasis:
        femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis,
                     1:DG_AffineLinearOnSimplexWithNodalBasis}
    else:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis,
                     1:DG_AffineP1_OnSimplexWithMonomialBasis}
elif spaceOrder == 2:
    if lagrangeBasis:
        femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis,
                      1:DG_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis,
                     1:DG_AffineP2_OnSimplexWithMonomialBasis}
elif spaceOrder == 3:
    assert not lagrangeBasis
    femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis,
                 1:DG_AffineP3_OnSimplexWithMonomialBasis}
else:
    femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
                 1:DG_AffineP0_OnSimplexWithMonomialBasis}


if spaceOrder == 0:
    elementQuadrature = SimplexLobattoQuadrature(nd,1)
    elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
else:
    elementQuadrature = SimplexGaussQuadrature(nd,2*spaceOrder+1)
    elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#if spaceOrder == 0:
#    nn=101

nn = 31#2**(8-spaceOrder)+1
nLevels = 1

subgridError = ShallowWater_CFL(coefficients,nd,g)

massLumping=False

shockCapturing = None

numericalFluxType = RusanovNumericalFlux#ShallowWaterHLL_1D#RusanovNumericalFlux#ShallowWater_1D

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
usingSSPRKNewton=False#True
multilevelNonlinearSolver=Newton
levelNonlinearSolver = Newton#SSPRKNewton#Newton
#fullNewtonFlag = False

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc#LU

    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = PETSc#LU#MGM#PETSc#
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    #parallelPartitioningType = MeshParallelPartitioningTypes.node
    parallelPartitioningType = MeshParallelPartitioningTypes.element
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU


linTolFac = 0.001

archiveFlag = ArchiveFlags.EVERY_USER_STEP
