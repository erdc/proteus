from pyadh import *
from pyadh.default_n import *
from la_gauss_1d_p import *

parallel = True

timeOrder =2
nStagesTime = timeOrder

DT=None
runCFL =0.1#3.0e-2

timeIntegration = LinearSSPRKintegration
stepController=Min_dt_RKcontroller
nDTout = 10

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 21
nLevels = 1

subgridError = None

massLumping = False

#numericalFluxType = Advection_DiagonalUpwind
numericalFluxType = RusanovLDG#Advection_DiagonalUpwind_Diffusion_LDG
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG

shockCapturing = None

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
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True
#maxNonlinearIts =1
#maxLineSearches =0
tolFac = 0.0

nl_atol_res = 1e-8

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
