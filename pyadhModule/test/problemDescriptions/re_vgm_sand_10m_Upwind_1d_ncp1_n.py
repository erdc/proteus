from pyadh import *
from pyadh.default_n import *
from re_vgm_sand_10m_Upwind_1d_p import *

#timeIntegration = BackwardEuler
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=51#3
nLevels = 1#5
DT = 1.0e-3
nDTout = 50#int(T/DT)

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-6

maxNonlinearIts = 100#0

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'point-eval'}
