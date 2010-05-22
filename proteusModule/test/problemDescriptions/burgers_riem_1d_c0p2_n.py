from pyadh import *
from pyadh.default_n import *
from burgers_riem_1d_p import *

#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF
stepController = FLCBDF_controller_sys
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3
runCFL=None#0.001

DT=None
nDTout = 1#int(T/DT)
print "nDTout",nDTout
femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd,5)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

#nn=3
#nLevels = 1
nn=26
nLevels=1
#subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)



massLumping=False
#massLumping=True

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
shockCapturing = ResGradDelayLag_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False,nStepsToDelay=5)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=True)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False)
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 25
maxLineSearches = 5

fullNewtonFlag = True

tolFac = 0.0#0.001

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None
