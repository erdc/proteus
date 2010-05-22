from pyadh import *
from pyadh.default_n import *
from twp_darcy_fcff_riem_1d_p import *

timeIntegrator = ForwardIntegrator
#timeIntegration = BackwardEuler
#stepController = FixedStep
#type of time integration formula
#timeIntegration = ForwardEuler_A
#general type of integration (Forward or to SteadyState)
stepController = FLCBDF_controller_sys
timeIntegration = FLCBDF
rtol_u[0] = 1.0e-3
rtol_u[1] = 1.0e-3
atol_u[0] = 1.0e-3
atol_u[1] = 1.0e-3
#runCFL = 1000.0
#runCFL = 20.0
runCFL=None

#DT=None
DT=1.0e1
nDTout = int(T/DT)
#nDTout=200
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,4)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
#nLevels = 1
nn=101
nLevels=1
#subgridError = None
#subgridError = FFDarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)


massLumping=False
#massLumping=True

shockCapturing = None
#shockCapturing = ResGradFFDarcy_SC(coefficients,nd,shockCapturingFactor=0.9,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=True)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=False)
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)

#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
#maxNonlinearIts = 25
maxNonlinearIts = 10
maxLineSearches = 10

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None
