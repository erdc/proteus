from pyadh import *
from pyadh.default_n import *
from kEpsilon_step_2d_p import *

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    runCFL = 0.1#None
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    atol_u[0] = 1.0e-2
    rtol_u[0] = 1.0e-2
    atol_u[1] = 1.0e-2
    rtol_u[1] = 1.0e-2

#runCFL = 0.5
#DT = None
#nDTout = 20

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,space_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,space_quad_order)

#nn=3
#triangleOptions = "Aq30Dena%f" % (0.15**2 / 6.0)
#triangleOptions = "Aq30Den"
#nLevels = 1


#subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

numericalFluxType = None

#shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
shockCapturing = ResGradDelayLag_SC(coefficients,nd,shockCapturingFactor=0.2,lag=False,nStepsToDelay=4)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
maxNonlinearIts =10

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True
