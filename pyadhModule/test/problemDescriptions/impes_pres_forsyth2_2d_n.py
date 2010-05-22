from pyadh import *
from pyadh.default_n import *
from impes_pres_forsyth2_2d_p import *

#type of time integration formula
timeIntegration = BackwardEuler
#timeIntegration = ForwardEuler_A
#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator

#runCFL = 1000.0
#runCFL = 20.0
runCFL=None

#DT=None
DT=1.0e1 
nDTout = int(T/DT)
#nDTout=200
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
#nLevels = 1
nn=3
nLevels=3
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)
subgridError = None


massLumping=False
#massLumping=True

#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.99,lag=True)
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)
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

conservativeFlux = {0:'point-eval'} #{0:'pwl'}
conservativeFlux = {0:'pwl-bdm'}
