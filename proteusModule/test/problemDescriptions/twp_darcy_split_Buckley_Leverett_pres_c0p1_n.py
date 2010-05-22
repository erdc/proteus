from pyadh import *
from pyadh.default_n import *
from twp_darcy_split_Buckley_Leverett_pres_p import *

#type of time integration formula
timeIntegration = NoIntegration

#runCFL = 1000.0
#runCFL = 20.0
runCFL=None

DT=None
#DT=5.0e-1
nDTout = 1#int(T/DT)
#nDTout=200
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
#nLevels = 1
nn=51
nLevels=1

subgridError = None
massLumping=False

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.99,lag=False)


multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
#maxNonlinearIts = 25

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

#conservativeFlux = {0:'pwl-bdm'}
#conservativeFlux = {0:'point-eval'}
conservativeFlux = {0:'pwl'}#'pwc'
