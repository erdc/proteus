from pyadh import *
from pyadh.default_n import *
from la_gauss_1d_p import *

#BackwardEuler
#timeIntegration = BackwardEuler_cfl
#stepController = Min_cfl_controller
DT = 1.0e-3
runCFL = 0.2
nDTout = 11
#timeIntegration = FLCBDF
#stepController = FLCBDF_controller
timeIntegration = VBDF
stepController = Min_dt_cfl_controller
timeOrder=2
#stepController = GustafssonFullNewton_dt_controller
#atol_u[0] = 1.0e-6
#rtol_u[0] = 1.0e-6
#timeOrder= 5

#now try simple error control?
#timeIntegrator  = testStuff.AdaptiveForwardIntegrator
#timeIntegration = testStuff.AdaptiveBackwardEuler

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

quad_order = 3
elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

nn = 51#101
nLevels = 1
#DT = 0.003125*0.5
#DT = None #use target cfl
#DT = 1.0e-6
#nDTout = 1

subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='1')


numericalFluxType = None
shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1)


multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

