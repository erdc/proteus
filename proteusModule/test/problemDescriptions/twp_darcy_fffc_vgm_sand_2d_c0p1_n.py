from pyadh import *
from pyadh.default_n import *
from twp_darcy_fffc_vgm_sand_2d_p import *


#timeIntegrator = ForwardIntegrator
#timeIntegration = BackwardEuler
#stepController = FixedStep
#DT=1.0
#runCFL=0.1
#runCFL=None
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-2
atol_u[0] = 1.0e-2
rtol_u[1] = 1.0e-2
atol_u[1] = 1.0e-2
#DT = None
nDTout = 50#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=81
#nLevels = 1

massLumping=False

subgridError = FFDarcyFC_ASGS(coefficients,nd,stabFlag='1',lag=True)
shockCapturing = ResGradFFDarcy_SC(coefficients,nd,shockCapturingFactor=0.25,lag=True)


multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts = 25
maxLineSearches = 10#0

fullNewtonFlag = True

tolFac = 0

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU
#pick number of layers to use in overlap 
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
nLayersOfOverlapForParallel = 1
#type of partition
parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element
numericalFluxType = DarcyFCFF_IIPG_exterior

linTolFac = 0.0001

conservativeFlux = {0:'pwl',1:'pwl'}
