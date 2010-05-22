from pyadh import *
from pyadh.default_n import *
from poisson_2c_ss_1d_p import *

timeIntegration = NoIntegration
nDTout = 1
parallel = True

femSpaces = dict((i,NC_AffineLinearOnSimplexWithNodalBasis) for i in range(coefficients.nc))


elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 21#5
nLevels = 1#5

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
    nLayersOfOverlapForParallel = 1
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None#{0:'p1-nc'} 
