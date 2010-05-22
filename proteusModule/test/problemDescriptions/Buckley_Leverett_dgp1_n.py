from pyadh import *
from pyadh.default_n import *
from Buckley_Leverett_p import *

timeOrder =2
nStagesTime = timeOrder

DT=None
runCFL = 0.3
limiterType =   TimeIntegration.DGlimiterP1Lagrange1d

timeIntegration = SSPRKPIintegration 
stepController=Min_dt_RKcontroller
nDTout = 1

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 26#41
nLevels = 1#1

subgridError = None

massLumping = False


numericalFluxType = RusanovNumericalFlux_Diagonal#

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
