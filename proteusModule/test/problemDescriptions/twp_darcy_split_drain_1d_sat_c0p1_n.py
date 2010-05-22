from pyadh import *
from pyadh.default_n import *
from twp_darcy_split_drain_1d_sat_p import *

#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator
#type of time integration formula
#timeIntegration = BackwardEuler
# stepController = FixedStep
#DT=1.0e1 
#stepController = Min_dt_controller
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[1] = 1.0e-5
rtol_u[2] = 1.0e-5
atol_u[1] = 1.0e-5
atol_u[2] = 1.0e-5
#runCFL = 1000.0
#runCFL = 200.0
runCFL=None
DT = None
nDTout = 50#int(T/DT)
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
#nLevels = 1
nn=101
nLevels=1
subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)


massLumping=False
#massLumping=True

#shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=True)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False)
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 20
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.0#0.001

atol = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = PETSc#PETSc LU

levelLinearSolver = PETSc#PETSc LU
#pick number of layers to use in overlap 
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
nLayersOfOverlapForParallel = 1
#type of partition
parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior


linTolFac = 0.0001

conservativeFlux = None
