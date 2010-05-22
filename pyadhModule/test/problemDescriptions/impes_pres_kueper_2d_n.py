from pyadh import *
from pyadh.default_n import *
from impes_pres_kueper_2d_p import *

#type of time integration formula
timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_A
#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator

#runCFL = 1000.0
#runCFL = 200.0
#runCFL=None
#runCFL=0.01


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
massLumping=False
shockCapturing = None

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
#maxNonlinearIts = 25

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None
conservativeFlux = {0:'point-eval'} #{0:'pwl'}
#conservativeFlux = {0:'pwl-bdm'}
