from proteus import *
from proteus.default_n import *
from la_gauss_2d_p import *

parallel = False

timeOrder = 2
nStagesTime = timeOrder
runCFL = 0.25
DT = None
nDTout = 10

usingSSPRKNewton = True
timeIntegration = LinearSSPRKintegration
stepController=Min_dt_RKcontroller
#limiterType = TimeIntegration.DGlimiterP1Lagrange2d#

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
#SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
#SimplexLobattoQuadrature(nd-1,1)

nn=21
nLevels = 1

subgridError = None

shockCapturing = None

numericalFluxType = Advection_DiagonalUpwind

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = SSPRKNewton#Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8
maxNonlinearIts =1

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = PETSc#LU

    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
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


linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
