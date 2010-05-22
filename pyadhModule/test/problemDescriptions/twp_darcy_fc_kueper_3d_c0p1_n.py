from pyadh import *
from pyadh.default_n import *
from twp_darcy_fc_kueper_3d_p import *

timeIntegrator = ForwardIntegrator
#type of time integration formula
#timeIntegration = BackwardEuler
#stepController = FixedStep
#DT = 1.0e1
#stepController = FixedStep
#stepController = FixedStep
#general type of integration (Forward or to SteadyState)
timeIntegration = FLCBDF_TwophaseDarcy_fc
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-3
rtol_u[1] = 1.0e-3
atol_u[0] = 1.0e-3
atol_u[1] = 1.0e-3
#runCFL=None
DT=None
nDTout = 50#int(T/DT)
#nDTout=200
print "nDTout",nDTout
print "T= ",T
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

triangleOptions = "pAfena%e" % ((top*right*top*0.1/6.0)*0.001,)
nLevels = 1

subgridError = None
#subgridError = DarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=True)


massLumping=False

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)#0.25 mostly
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)

#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
#maxNonlinearIts = 25
maxNonlinearIts = 20
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc
#pick number of layers to use in overlap 
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
nLayersOfOverlapForParallel = 1
#type of partition
parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element
numericalFluxType = DarcyFC_IIPG_exterior

linTolFac = 0.0001

conservativeFlux = None
