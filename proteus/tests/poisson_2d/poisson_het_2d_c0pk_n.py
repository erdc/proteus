from __future__ import absolute_import
from builtins import range
from proteus import *
from proteus.default_n import *
try:
    from .poisson_het_2d_p import *
except:
    from poisson_het_2d_p import *

parallel = False
polynomial_order = 2
timeIntegration = NoIntegration
nDTout = 1

if polynomial_order == 2:
    femSpaces = dict((i,C0_AffineQuadraticOnSimplexWithNodalBasis) for i in range(nc))
else:    
    femSpaces = dict((i,C0_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 21
nLevels = 2
if parallel:
    nn = (nn-1)*2**(nLevels-1)+1
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

class Exterior_StrongFlux(DoNothing):
    useStrongDirichletConstraints=True

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
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    linearSolverConvergenceTest= 'r'

linTolFac = 0.0
l_atol_res = 1.0e-10

if polynomial_order == 2:
    cfluxtag = 'pwl-bdm'
else:
    cfluxtag  = 'pwl'#'pwl-bdm'#'sun-rt0','sun-gs-rt0','pwc','pwl','pwl-bdm','point-eval'
conservativeFlux =  dict((i,cfluxtag) for i in range(nc))
#need this for sun-wheeler-gs
if cfluxtag == 'sun-gs-rt0':
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG_exterior

