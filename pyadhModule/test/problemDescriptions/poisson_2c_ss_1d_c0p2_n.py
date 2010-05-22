from pyadh import *
from pyadh.default_n import *
from poisson_2c_ss_1d_p import *

parallel = True
timeIntegration = NoIntegration
nDTout = 1


femSpaces = dict((i,C0_AffineQuadraticOnSimplexWithNodalBasis) for i in range(coefficients.nc))


elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 21#2**2+1
nLevels = 1#5

subgridError = None

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc

    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = PETSc
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU


linearSmoother = StarILU#GaussSeidel#Jacobi#StarILU

linTolFac = 1.0e-6

multigridCycles = 3

preSmooths = 3

postSmooths = 3

cfluxtag  = 'point-eval'
conservativeFlux =  dict((i,cfluxtag) for i in range(coefficients.nc))

