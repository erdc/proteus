from pyadh import *
from pyadh.default_n import *
from Buckley_Leverett_Liu_2d_p import *

useBackwardEuler = False
timeIntegrator = ForwardIntegrator
if useBackwardEuler:
   timeIntegration = BackwardEuler_cfl
   stepController = Min_dt_controller
   runCFL = 0.1
   DT = None 
else:
    #timeIntegration = FLCBDF
    #stepController = FLCBDF_controller_sys
    timeIntegration = VBDF
    stepController = GustafssonFullNewton_dt_controller
    rtol_u[0] = 1.0e-2
    rtol_u[0] = 1.0e-2
    runCFL = None


nDTout = 50

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 1
triangleOptions = "q30Dena0.005A"
nLevels = 1#3



subgridError = None
#subgridError = Advection_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='1',lag=True)
subgridError = AdvectionDiffusionReactionTransientSubscales_ASGS(coefficients,nd,stabFlag='1',lag=True,trackSubScales=True,useHarariDirectly=False,
                                                                 limit_tau_t=True,tau_t_limit_max=0.9)

massLumping = False

numericalFluxType = None
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)

ultilevelNonlinearSolver  = Newton#NLNI

levelNonlinearSolver = Newton#

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-7

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

