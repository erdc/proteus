from __future__ import absolute_import
from builtins import range
from proteus import *
from proteus.default_n import *
try:
    from .poisson_het_2d_p import *
except:
    from poisson_het_2d_p import *

parallel = False
numerical_flux_flag = 'SIPG'
polynomial_order = 2


timeIntegration = NoIntegration
nDTout = 1

if polynomial_order == 2:
    femSpaces = dict((i,DG_AffineQuadraticOnSimplexWithNodalBasis) for i in range(nc))
else:
    femSpaces = dict((i,DG_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 11
nLevels = 1
if parallel:
    nLevels = 1
    nn = nn*2**(nLevels-1)
    nLevels = 1

subgridError = None

shockCapturing = None

if numerical_flux_flag == 'SIPG':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG
elif numerical_flux_flag == 'IIPG':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
elif numerical_flux_flag == 'LDG':
    numericalFluxType = Diffusion_LDG
else:
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_NIPG

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
    #matrix = np.array
linTolFac = 0.0

if polynomial_order == 2:
   cfluxtag  = 'dg-point-eval' #'dg-point-eval','dg'
else:
   cfluxtag  = 'dg'
conservativeFlux = dict((i,cfluxtag) for i in range(nc))
