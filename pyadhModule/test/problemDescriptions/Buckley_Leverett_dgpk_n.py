from pyadh import *
from pyadh.default_n import *
from Buckley_Leverett_p import *

parallel = True
timeOrder =3
nStagesTime = timeOrder

DT=None
runCFL = 0.125 #0.1
limiterType = TimeIntegration.DGlimiterPkMonomial1d

#note about initial condition
#need something smoothed out right now for dgpk because default projection of step data
#causes oscillations


#mwf BackwardEuler 
timeIntegration = SSPRKPIintegration 
stepController=Min_dt_RKcontroller
nDTout = 1

#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
#highest can go with current quadrature
femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 41
nLevels = 1

subgridError = None

massLumping = False

numericalFluxType = RusanovNumericalFlux_Diagonal#testStuff.RusanovNumericalFlux_Diagonal_Diffusion_LDG#

shockCapturing = None

multilevelNonlinearSolver  = NLNI
usingSSPRKNewton = True
levelNonlinearSolver = SSPRKNewton#Newton

nonlinearSmoother = NLGaussSeidel

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
