from proteus import *
from proteus.default_n import *
from .poisson_het_2d_p import *

parallel = False

timeIntegration = NoIntegration
nDTout = 1

femSpaces = dict((i,NC_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 21
nLevels = 2
if parallel:
    nLevels = 1

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
maxNonlinearIts = 1

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py
    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = KSP_petsc4py
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #to allow multiple models to set different ksp options
    #linear_solver_options_prefix = 'poisson_'
    linearSmoother = None
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    linearSolverConvergenceTest= 'r'


linTolFac = 0.0
l_atol_res = 1.0e-10

cfluxtag  = 'p1-nc'
conservativeFlux =  dict((i,cfluxtag) for i in range(nc))
