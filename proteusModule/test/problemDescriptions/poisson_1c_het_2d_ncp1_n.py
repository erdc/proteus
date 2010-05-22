from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_2d_p import *

timeIntegration = NoIntegration
nDTout = 1
parallel = True

femSpaces = dict((i,NC_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))
elementQuadrature =SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 11#5
nLevels = 3#5

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton# NLNI

levelNonlinearSolver = Newton

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
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
else:
    multilevelLinearSolver = LU#NI#MGM

    levelLinearSolver = LU#MGM

linearSmoother = GaussSeidel

linTolFac = 0.001

cfluxtag  = 'p1-nc'
conservativeFlux =  dict((i,cfluxtag) for i in range(nc))

