from pyadh import *
from pyadh.default_n import *
from la_cone_2d_p import *

timeOrder = 2
nStagesTime = timeOrder
runCFL = 0.3
timeIntegration = SSPRKPIintegration
limiterType = TimeIntegration.DGlimiterDurlofskyP1Lagrange2d

stepController=Min_dt_RKcontroller
DT = None
nDTout = 20

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)#SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)#SimplexLobattoQuadrature(nd-1,1)

nn = 41
nLevels = 1#4

subgridError = None

massLumping = False

shockCapturing = None

numericalFluxType = Advection_DiagonalUpwind

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = SSPRKNewton#Newton

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
