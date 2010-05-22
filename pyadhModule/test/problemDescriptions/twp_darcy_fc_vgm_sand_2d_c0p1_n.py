from pyadh import *
from pyadh.default_n import *
from twp_darcy_fc_vgm_sand_2d_p import *

#type of time integration formula
#timeIntegration = BackwardEuler
#stepController = FixedStep
#DT=1.0e2
#nDTout = int(T/DT)
#timeIntegration = OuterTheta
#timeIntegration = ForwardEuler_A
#general type of integration (Forward or to SteadyState)
#timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF_TwophaseDarcy_fc
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-2
atol_u[0] = 1.0e-2
rtol_u[1] = 1.0e-2
atol_u[1] = 1.0e-2
DT = None
nDTout = 50#int(T/DT)
runCFL=None

print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
#nLevels = 1
nn=51
nLevels=1
subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
subgridError = DarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=True)


massLumping=False
#massLumping=True

#shockCapturing = None
#shockCapturing = ResGradFFDarcy_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)
shockCapturing = ResGradDelayLag_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False,nStepsToDelay=5)#0.25 mostly
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)

multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 10
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.1

nl_atol_res = 1.0e-7#1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = PETSc#LU

levelLinearSolver = PETSc#LU
#pick number of layers to use in overlap 
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
nLayersOfOverlapForParallel = 1
#type of partition
parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element
numericalFluxType = DarcyFC_IIPG_exterior


linTolFac = 0.0001

#conservativeFlux = {0:'point-eval'} #{0:'pwl'}
conservativeFlux = None
