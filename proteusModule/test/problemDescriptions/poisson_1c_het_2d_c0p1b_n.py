from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_2d_p import *

parallel = False
timeIntegration = NoIntegration
nDTout = 1

dgBubbles = True#False
if dgBubbles:
    femSpaces = dict((i,C0_AffineP1P0BubbleOnSimplexWithNodalBasis) for i in range(nc))
else:
    femSpaces = dict((i,C0_AffineP1BubbleOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 41#21
nLevels = 1#7

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton#NLNI
levelNonlinearSolver = Newton
maxNonlinearIts = 1

fullNewtonFlag = True

tolFac = 1.0e-8

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
    nLayersOfOverlapForParallel = 2
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU#MGM#
if dgBubbles:
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
linearSmoother = StarILU#GaussSeidel#Jacobi#StarILU

linTolFac = 1.0e-6

multigridCycles = 3

cfluxtag  = 'dg-point-eval'#'sun-rt0','sun-gs-rt0','pwc','pwl','pwl-bdm','point-eval'
conservativeFlux =  dict((i,cfluxtag) for i in range(nc))

preSmooths = 3

postSmooths = 3


