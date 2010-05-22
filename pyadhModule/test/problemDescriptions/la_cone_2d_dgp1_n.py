from pyadh import *
from pyadh.default_n import *
from la_cone_2d_p import *

timeOrder = 3
nStagesTime = timeOrder
runCFL = 0.1
timeIntegration = SSPRKPIintegration 
stepController=Min_dt_RKcontroller

DT = None
nDTout = 20

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)#SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)#SimplexLobattoQuadrature(nd-1,1)

nn = 17
nLevels = 1

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
