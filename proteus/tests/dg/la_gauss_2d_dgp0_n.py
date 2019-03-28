from proteus import *
from proteus.default_n import *
from la_gauss_2d_p import *

timeOrder = 1
nStagesTime = timeOrder
runCFL = 0.45
DT = None
nDTout = 10

#BackwardEuler, ForwardEuler
timeIntegration = LinearSSPRKintegration
stepController=Min_dt_RKcontroller

femSpaces = {0:DG_Constants}

# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nLevels = 4

subgridError = None

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
