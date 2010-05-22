from pyadh import *
from pyadh.default_n import *
from twp_darcy_fc_hole_2x4m_2d_p import *

timeIntegrator = ForwardIntegrator
useFLCBDF = True
useVBDF = False

if useFLCBDF:
    #timeIntegration = FLCBDF_TwophaseDarcy_fc
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller
    runCFL=None
    DT=None
    nDTout = 50#int(T/DT)
elif useVBDF:
    timeIntegration = VBDF
    timeOrder =2
    stepController = GustafssonFullNewton_dt_controller
    nonlinearConvergenceRate_ref = 0.8
    useInitialGuessPredictor = False#True
    stepExact = True
    runCFL=None
    DT=None
    nDTout = 50#int(T/DT)
else:
    timeIntegration = BackwardEuler_cfl
    stepController =  Min_dt_controller
    DT = None
    nDTout = 50
    runCFL = 0.1

#type of time integration formula
#timeIntegration = BackwardEuler_cfl
#stepController = Min_dt_controller
#runCFL = 0.3
#timeIntegration = FLCBDF
#stepController  = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
#general type of integration (Forward or to SteadyState)
#timeIntegration = FLCBDF_TwophaseDarcy_fc
#stepController = FLCBDF_controller

#timeIntegration = VBDF
#timeOrder =1
#stepController = GustafssonFullNewton_dt_controller
#nonlinearConvergenceRate_ref = 0.8
#useInitialGuessPredictor = True
#stepExact = True
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

rtol_u[0] = 1.0e-3
rtol_u[1] = 1.0e-3
atol_u[0] = 1.0e-3
atol_u[1] = 1.0e-3
#runCFL=None
DT=None
nDTout = 50#int(T/DT)

print "nDTout",nDTout
print "T= ",T
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
triangleOptions = "q30Dena0.005A"   
#triangleOptions = "q30Dena0.005A"   
#triangleOptions = "pq30Dena0.001A"
#triangleOptions = "pq30Dena0.0025A"
nLevels = 1

subgridError = None
subgridError = DarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=False)


massLumping=False

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)#0.25 mostly
#shockCapturing = ResGradDelayLag_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False,nStepsToDelay=5)#0.25 mostly

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts = 10
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-5

matrix = SparseMatrix

multilevelLinearSolver = LU#PETSc#PETSc LU

levelLinearSolver = LU#PETSc#PETSc LU
#pick number of layers to use in overlap 
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
#nLayersOfOverlapForParallel = 2
#type of partition
#parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element
#numericalFluxType = DarcyFC_IIPG_exterior
#numericalFluxType = StrongDirichletFactory(fluxBoundaryConditions)#None
#numericalFluxType = StrongDirichletFactory(fluxBoundaryConditions)#None

linTolFac = 0.0001

conservativeFlux = None
