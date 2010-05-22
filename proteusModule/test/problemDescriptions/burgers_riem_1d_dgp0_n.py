from pyadh import *
from pyadh.default_n import *
from burgers_riem_1d_p import *

timeOrder =2
nStagesTime = timeOrder

DT=None
runCFL = 0.3 #0.1 0.8

timeIntegration = SSPRKPIintegration 
stepController=Min_dt_RKcontroller
nDTout = 10

femSpaces = {0:DG_Constants}

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 101
nLevels = 1

subgridError = None

massLumping = False

numericalFluxType = Advection_DiagonalUpwind

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
