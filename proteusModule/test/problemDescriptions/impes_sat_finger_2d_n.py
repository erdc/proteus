from pyadh import *
from pyadh.default_n import *
from impes_sat_finger_2d_p import *

#type of time integration formula
timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_A
#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator

#runCFL = 1000.0
#runCFL = 200.0
runCFL=None
#runCFL=0.1

DT=1.0e1 
#DT=None
nDTout = int(T/DT)
#nDTout=100
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
#nLevels = 1
nn=3
nLevels=6
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = None


massLumping=False
#massLumping=True

shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.99,lag=False)
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
#shockCapturing = None

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 25

fullNewtonFlag = True

tolFac = 0.001

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None
