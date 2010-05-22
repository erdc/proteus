from pyadh import *
from pyadh.default_n import *
from re_het_column_8_3d_p import *

useFLCBDF = False
useGustafsson = True#False#True
atol_u[0] = 1.0e-2#1.0e-6
rtol_u[0] = 1.0e-2#1.0e-6
tnList=[0.0,T*1.0e-4,T]
nDTout = len(tnList)

if useFLCBDF:
    timeIntegration = FLCBDF
    stepController  = FLCBDF_controller
    systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
elif useGustafsson:
    #timeIntegration = BackwardEuler
    timeIntegration = VBDF
    timeOrder = 2
    stepController = GustafssonFullNewton_dt_controller#FixedStep
    #for controlling time stepping
    nonlinearConvergenceRate_ref = 0.4#0.6#0.3#0.3#0.2
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
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#elementQuadrature = SimplexGaussQuadrature(nd,4)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nLevels = 1
triangleOptions="VpAq1.25ena%f" % (((0.1*L[2])**3)/6.0,)

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

shockCapturing = None
#shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True,nStepsToDelay=1)

numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
#multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
multilevelNonlinearSolver = Newton

#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

#nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-5

maxNonlinearIts = 25
maxLineSearches =25

matrix = SparseMatrix
#matrix = numpy.array

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = PETSc
#levelLinearSolver = MGM

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

#conservativeFlux = {0:'pwl'}
#conservativeFlux = {0:'point-eval'}
