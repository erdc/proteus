from pyadh import *
from pyadh.default_n import *
from la_hauke_sangalli_1d_p import *

timeOrder = 3
nStagesTime = timeOrder
timeIntegration = LinearSSPRKintegration
stepController  = Min_dt_RKcontroller
runCFL = 0.1
 
nLevels = 1
nn = haukeSangalli_nn
nDTout = int(T/haukeSangalli_dt)

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


subgridError = None
numericalFluxType = Advection_DiagonalUpwind#_Diffusion_IIPG

multilevelNonlinearSolver  = NLNI

maxNonlinearIts =1

maxLineSearches =0

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True


archiveFlag = ArchiveFlags.EVERY_USER_STEP
