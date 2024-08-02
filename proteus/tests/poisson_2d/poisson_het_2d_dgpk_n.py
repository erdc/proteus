from proteus import *
from proteus.default_n import *
from proteus import defaults
defaults.reset_default_n()
try:
    from .poisson_het_2d_p import *
except:
    from poisson_het_2d_p import *

parallel = True
direct=False
polynomial_order = 2
numerical_flux_flag = 'SIPG'

timeIntegration = NoIntegration
nDTout = 1

if polynomial_order == 2:
    femSpaces = dict((i,DG_AffineQuadraticOnSimplexWithNodalBasis) for i in range(nc))
else:
    femSpaces = dict((i,DG_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 81
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
maxNonlinearIts = 1

fullNewtonFlag = True

tolFac = 0.0
nl_atol_res = 1.0e-8
linTolFac = 0.0
l_atol_res = 1.0e-9

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    #pick number of layers to use in overlap and type of partition
    if numerical_flux_flag == 'LDG':
        nLayersOfOverlapForParallel = 2
    else:
        nLayersOfOverlapForParallel = 1
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.element
    linearSolverConvergenceTest= 'r-true'
    from petsc4py import PETSc
    OptDB = PETSc.Options()
    OptDB.clear()
    if direct:
        OptDB.setValue('ksp_type','preonly')
        OptDB.setValue('pc_type','lu')
        OptDB.setValue('pc_factor_type','superlu_dist')
    else:
        OptDB.setValue('ksp_type','cg')
        OptDB.setValue('pc_asm_type','basic')
        OptDB.setValue('pc_asm_overlap',2)
        OptDB.setValue('sub_ksp_type','preonly')
        OptDB.setValue('sub_pc_type','lu')
        OptDB.setValue('sub_pc_factor_type','superlu')
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    linearSolverConvergenceTest= 'r'

conservativeFlux = None

hex = False
quad = False
