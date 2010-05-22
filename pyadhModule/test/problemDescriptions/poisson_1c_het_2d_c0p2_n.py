from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_2d_p import *

parallel = False
timeIntegration = NoIntegration
nDTout = 1


femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd,6)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,6)

nnx = 11#2**2+1
nny = 11
nLevels = 3#5

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
    multilevelLinearSolver = PETSc#LU

    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = PETSc#LU#MGM#PETSc#
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    #parallelPartitioningType = MeshParallelPartitioningTypes.node
    parallelPartitioningType = MeshParallelPartitioningTypes.element
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    #mwf debug
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior


linearSmoother = StarILU#GaussSeidel#Jacobi#StarILU

linTolFac = 1.0e-6

multigridCycles = 3

preSmooths = 3

postSmooths = 3

# #if using sun-gs-rt0 need weak dirichlet conditions
# numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
conservativeFlux =  {0:'pwl-bdm'}#{0:'sun-rt0'}#{0:'pwc'}#{0:'pwl'}#{0:'pwl-bdm'}#{0:'point-eval'}#

