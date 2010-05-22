from pyadh import *
from pyadh.default_n import *
from sw_Sivakumaran_flume_2d_p import *

parallel = True

#runCFL=0.33
runCFL=0.15
timeOrder = 2
nStagesTime=timeOrder
timeIntegration = SSPRKPIintegration
nDTout = 101
tnList = [0.0 + i*T/(float(nDTout)-1.) for i in range(nDTout)]
limiterType = DGlimiterDurlofskyP1Lagrange2d_Sw#DGlimiterP1Lagrange2d#DGlimiterDurlofskyP1Lagrange2d

stepController=Min_dt_RKcontroller
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

#nDTout=1001
femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis,
             1:DG_AffineLinearOnSimplexWithNodalBasis,
             2:DG_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elementQuadrature = SimplexGaussQuadrature(nd,2)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

nnx=101
nny=21
nLevels = 1

subgridError = ShallowWater_CFL(coefficients,nd,g)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

massLumping=False

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#numericalFluxType = Advection_DiagonalUpwind
numericalFluxType = ShallowWater_2D

multilevelNonlinearSolver  = NLNI

#usingSSPRKNewton=True
#levelNonlinearSolver = SSPRKNewton#Newton
#fullNewtonFlag = False
usingSSPRKNewton=False
levelNonlinearSolver = Newton
tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc#LU
    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = PETSc#LU#MGM#PETSc#
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.element
    #parallelPartitioningType = MeshParallelPartitioningTypes.node

else:
    multilevelLinearSolver = LU#NI#MGM

    levelLinearSolver = LU#MGM


linTolFac = 0.001

#conservativeFlux = {0:'pwl'}

archiveFlag = ArchiveFlags.EVERY_USER_STEP
