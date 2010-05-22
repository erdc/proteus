from pyadh import *
from pyadh.default_n import *
from twp_darcy_fc_kueper_2d_p import *

timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator
runCFL=None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
massLumping=True
shockCapturing = None

#subgridError = FFDarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=False)
#massLumping=False
#shockCapturing = ResGradFFDarcy_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 25

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None
