from proteus import *
from proteus.default_n import *
from poisson_het_2d_p import *

parallel = True
numerical_flux_flag = 'SIPG'


timeIntegration = NoIntegration
nDTout = 1

femSpaces = dict((i,DG_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 11
if parallel:
    nLevels = 1
else:
    nLevels = 4

subgridError = None

shockCapturing = None

if numerical_flux_flag == 'SIPG':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG
elif numerical_flux_flag == 'NIPG':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_NIPG
elif numerical_flux_flag == 'LDG':
    numericalFluxType = Diffusion_LDG
else:
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

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
    #pick number of layers to use in overlap and type of partition
    if numerical_flux_flag == 'LDG':
        nLayersOfOverlapForParallel = 2
    else:
        nLayersOfOverlapForParallel = 1
    parallelPartitioningType = MeshParallelPartitioningTypes.element
    
    #to allow multiple models to set different ksp options
    #linear_solver_options_prefix = 'poisson_'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    linearSolverConvergenceTest= 'r'

linTolFac = 0.0

cfluxtag  = 'dg' #'dg-point-eval','dg'
conservativeFlux = dict((i,cfluxtag) for i in range(nc))
