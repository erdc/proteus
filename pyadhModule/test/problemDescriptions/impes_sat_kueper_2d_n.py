from pyadh import *
from pyadh.default_n import *
from impes_sat_kueper_2d_p import *

#type of time integration formula
timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_A
#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator

#runCFL = 1000.0
#runCFL = 200.0
#runCFL=None
#runCFL=0.01

#nDTout=1
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)


massLumping=False
#massLumping=True

shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 25

fullNewtonFlag = True

tolFac = 0.1

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None
