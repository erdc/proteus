from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_3d_p import *
parallel = True
timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 5#11
nLevels = 1#5

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-4

nl_atol_res = 1.0e-4

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc#NI#LU

    levelLinearSolver = PETSc#MGM#LU

    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    #parallelPartitioningType = MeshParallelPartitioningTypes.node
    parallelPartitioningType = MeshParallelPartitioningTypes.element
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU
    
linTolFac = 0.001

#conservativeFlux = {0:'pwl-bdm'} #{0:'pwl-bdm'} #{0:'point-eval'}#{0:'p1-nc'}#{0:'pwc'} #



