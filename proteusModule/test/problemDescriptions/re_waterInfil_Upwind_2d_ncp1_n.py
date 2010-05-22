from pyadh import *
from pyadh.default_n import *
from re_waterInfil_Upwind_2d_p import *

#timeIntegration = BackwardEuler
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3
DT = None
femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
#triangleOptions += "A"
triangleOptions = "q30Dena0.005A"   
#triangleOptions = "pq30Dena0.001A"
#triangleOptions = "pq30Dena0.0025A"
nLevels = 1
nDTout = 1#int(T/DT)

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)
#seems to have more under/overshoot 
#shockCapturing = GenericJaffre_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True,betaPower=0.1)
#delay lag not working well right now
#shockCapturing = ResGradDelayLag_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False,nStepsToDelay=5)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

nLayersOfOverlapForParallel = 2
parallelPartitioningType = MeshParallelPartitioningTypes.node

tolFac = 0.0

nl_atol_res = 1.0e-6

maxNonlinearIts = 20#0
maxLineSearches = 5

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'p1-nc'}
