from pyadh import *
from pyadh.default_n import *
from re_vgm_clay_30cm_1d_p import *


timeIntegrator = ForwardIntegrator
useFLCBDF = True#False#True
useGustafsson = False#True
atol_u[0] = 1.0e-2#1.0e-2#
rtol_u[0] = 1.0e-2#1.0e-2#
DT = 1.0e-5#None#0.025#1.0e-1/timeScale
tnList = [0.0,DT,T]; nDTout = 3
if useFLCBDF:
    timeIntegration = FLCBDF
    stepController  = FLCBDF_controller
    systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
elif useGustafsson:
    #timeIntegration = BackwardEuler
    timeIntegration = VBDF
    timeOrder = 1
    stepController = GustafssonFullNewton_dt_controller#FixedStep
    #for controlling time stepping
    nonlinearConvergenceRate_ref = 0.2#0.6#0.3#0.3#0.2
    useInitialGuessPredictor= True
    stepExact = True
    systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
else:
    timeIntegration = BackwardEuler
    stepController = HeuristicNL_dt_controller#FixedStep
    #nDTout = 1#int(T/DT)#int(T/DT) #100#int(T/DT)
    #for controlling time stepping
    nonlinearIterationsFloor = 4
    nonlinearIterationsCeil  = 8
    dtNLgrowFactor = 2
    dtNLreduceFactor = 0.5
    dtNLfailureReduceFactor = 0.5
    useInitialGuessPredictor= True
    stepExact = True

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexLobattoQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=51#3
nLevels = 1#5

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)
subgridError = AdvectionDiffusionReactionTransientSubscales_ASGS(coefficients,nd,stabFlag='2',lag=True,trackSubScales=True,useHarariDirectly=True)

massLumping = False

numericalFluxType = None

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-6

maxNonlinearIts = 10#0

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
