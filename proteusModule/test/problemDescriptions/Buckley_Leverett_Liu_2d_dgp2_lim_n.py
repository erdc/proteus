from pyadh import *
from pyadh.default_n import *
from Buckley_Leverett_Liu_2d_p import *

parallel = True#True

timeOrder =3
nStagesTime = timeOrder

runCFL = 0.1 #0.1

limiterType = TimeIntegration.DGlimiterP2Lagrange2d#

timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller
DT=None
nDTout = 1

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 1
nLevels = 1#3
triangleOptions="q30a0.005Den"
subgridError = None

massLumping = False

numericalFluxType =  RusanovNumericalFlux_Diagonal#

shockCapturing = None

multilevelNonlinearSolver  = Newton#NLNI

usingSSPRKNewton=True
levelNonlinearSolver = SSPRKNewton#
#usingSSPRKNewton=False
#levelNonlinearSolver = Newton#

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

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

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None


archiveFlag = ArchiveFlags.EVERY_USER_STEP
