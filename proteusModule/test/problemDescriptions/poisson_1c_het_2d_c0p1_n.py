from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_2d_p import *

parallel = True
timeIntegration = NoIntegration
nDTout = 1


femSpaces = dict((i,C0_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 21#3#41#21
nLevels = 1#5#7

    
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
    multilevelLinearSolver = KSP_petsc4py#PETSc#LU
    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = KSP_petsc4py#PETSc#LU#MGM#PETSc#
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    #parallelPartitioningType = MeshParallelPartitioningTypes.node
    parallelPartitioningType = MeshParallelPartitioningTypes.element
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
    #for true residual test, doesn't work for petsc4py
    #linearSolverConvergenceTest = 'r-true'
    #to allow multiple models to set different ksp options
    linear_solver_options_prefix = 'poisson_'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU#MGM#
    #numericalFluxType = None
    #numericalFluxType = StrongDirichletFactory(fluxBoundaryConditions)
    linearSolverConvergenceTest= 'r'#r-true'#'r'

linearSmoother = StarILU#GaussSeidel#Jacobi#StarILU

linTolFac = 0.0
l_atol_res = 1.0e-10

multigridCycles = 3

cfluxtag  = 'pwl'#'pwl-bdm'#'sun-rt0','sun-gs-rt0','pwc','pwl','pwl-bdm','point-eval'
conservativeFlux =  {0:'pwl',1:'pwl-ib-fix-0'}#dict((i,cfluxtag) for i in range(nc))
#need this for sun-wheeler-gs
if cfluxtag == 'sun-gs-rt0':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
preSmooths = 3

postSmooths = 3

