from pyadh import *
from pyadh.default_n import *
from twp_darcy_split_hole_2x4m_2d_pres_p import *

parallel = True
#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator
timeIntegration = NoIntegration

runCFL=None
DT=None
nDTout=50
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=3
nLevels = 1
triangleOptions = "q30Dena0.005A"
nLevels = 1

subgridError = None
massLumping=False

shockCapturing = None


#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts = 20
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix
if parallel:

    multilevelLinearSolver = PETSc#PETSc LU

    levelLinearSolver = PETSc#PETSc LU
    #pick number of layers to use in overlap 
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    nLayersOfOverlapForParallel = 2
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    numericalFluxType = Diffusion_IIPG_exterior
else:

    multilevelLinearSolver = LU

    levelLinearSolver = LU
     


linTolFac = 0.0001

#conservativeFlux = {0:'pwl-bdm'}
#conservativeFlux = {0:'point-eval'}
conservativeFlux = {0:'pwl'}
