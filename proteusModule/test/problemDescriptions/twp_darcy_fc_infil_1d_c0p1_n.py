from pyadh import *
from pyadh.default_n import *
from twp_darcy_fc_infil_1d_p import *

timeIntegrator = ForwardIntegrator
useFLCBDF = False
useVBDF = True
#type of time integration formula
#works with compressibility
#timeIntegration = BackwardEuler
#stepController = FixedStep
#DT = 1.0e-1
rtol_u[0] = 1.0e-4
rtol_u[1] = 1.0e-4
atol_u[0] = 1.0e-4
atol_u[1] = 1.0e-4
if useFLCBDF:
    timeIntegration = FLCBDF_TwophaseDarcy_fc
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller
    runCFL=None
    DT=None
    nDTout = 50#int(T/DT)
elif useVBDF:
    timeIntegration = VBDF
    timeOrder =2
    stepController = GustafssonFullNewton_dt_controller
    nonlinearConvergenceRate_ref = 0.2
    useInitialGuessPredictor = False#True
    stepExact = False#True
    runCFL=None
    DT=None
    nDTout = 50#int(T/DT)
else:
    timeIntegration = BackwardEuler_cfl
    stepController =  Min_dt_controller
    DT = 0.1#None
    nDTout = 50
    runCFL = 0.05
#nDTout=200
tnList=[0.,0.1]#
tnList.extend([i*T/nDTout for i in range(1,nDTout+1)])
print "nDTout",nDTout
print "T= ",T
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=101
nLevels=1
subgridError = None
#subgridError = DarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=True)


massLumping=False
#massLumping=True

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)#0.25 mostly

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts = 20
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

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
#numericalFluxType = DarcyFC_IIPG_exterior

linTolFac = 0.0001

conservativeFlux = None
