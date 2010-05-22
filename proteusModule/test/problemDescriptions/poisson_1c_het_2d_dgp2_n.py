from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_2d_p import *

parallel = False#True

timeIntegration = NoIntegration
nDTout = 1


femSpaces = dict((i,DG_AffineQuadraticOnSimplexWithNodalBasis) for i in range(nc))


elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 5
nLevels = 3
if parallel:
    nn = nn*2**(nLevels-1)
    nLevels = 1
subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5e2,lag=True)
#shockCapturing = ResGradJuanes_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)

numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG#Diffusion_LDG#

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-6

nl_atol_res = 1.0e-12

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



linearSmoother = StarILU#Jacobi#StarILU

linTolFac = 1.0e-8

cfluxtag  = 'dg-bdm' #'dg-point-eval','dg','dg-bdm'
conservativeFlux = dict((i,cfluxtag) for i in range(nc))

multigridCycles = 3

preSmooths = 3

postSmooths = 3

archiveFlag = ArchiveFlags.EVERY_USER_STEP
