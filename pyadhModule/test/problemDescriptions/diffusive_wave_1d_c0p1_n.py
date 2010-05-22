from pyadh import *
from pyadh.default_n import *
from diffusive_wave_1d_p import *


timeIntegrator = ForwardIntegrator
useFLCBDF = False    #True#False
useGustafsson = False#False#True#
atol_u[0] = 1.0e-4
rtol_u[0] = 1.0e-4
DT = 1.0e-3#None#0.025#1.0e-1/timeScale
tnList = [0.0,DT,T]; nDTout = 3
if useFLCBDF:
    timeIntegration = FLCBDF
    stepController  = FLCBDF_controller
    systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
elif useGustafsson:
    timeIntegration = VBDF
    timeOrder = 2
    stepController = GustafssonFullNewton_dt_controller
    #for controlling time stepping
    nonlinearConvergenceRate_ref = 0.2
    useInitialGuessPredictor= True
    stepExact = True
    systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
else:
    timeIntegration = BackwardEuler
    stepController = HeuristicNL_dt_controller#FixedStep
    #for controlling time stepping
    nonlinearIterationsFloor = 4
    nonlinearIterationsCeil  = 8
    dtNLgrowFactor = 2
    dtNLreduceFactor = 0.5
    dtNLfailureReduceFactor = 0.5
    useInitialGuessPredictor= True
    stepExact = True

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elementQuadrature = SimplexGaussQuadrature(nd,3)
#
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=101
nLevels = 1

subgridError = None

masslumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLStarILU

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6#1.0e-8

maxNonlinearIts = 10#1001
maxLineSearches =5

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

#conservativeFlux = {0:'pwl'}
