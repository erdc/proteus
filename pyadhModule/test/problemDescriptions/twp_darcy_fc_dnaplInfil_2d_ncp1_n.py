from pyadh import *
from pyadh.default_n import *
from twp_darcy_fc_dnaplInfil_2d_p import *

timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator
#runCFL=0.1
DT = 1.0e2#1.0e-2
#timeIntegration = FLCBDF
#stepController = FLCBDF_controller
#rtol_u[0] = 1.0e-3
#atol_u[0] = 1.0e-3
#rtol_u[1] = 1.0e-3
#atol_u[1] = 1.0e-3
#DT = None
nDTout = 50#int(T/DT)

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis,
             1:NC_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=41
nLevels = 1



massLumping=False
subgridError = None
shockCapturing = None


multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 10#025
maxLineSearches = 10#0

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = {0:'p1-nc',1:'p1-nc'}
