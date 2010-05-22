from pyadh import *
from pyadh.default_n import *
from la_shock_1d_p import *

#mwf add
timeOrder =3
nStagesTime = timeOrder
runCFL = 0.9
DT = None
nDTout = 20
timeIntegration = SSPRKPIintegration
limiterType = TimeIntegration.DGlimiterP1Lagrange1d

stepController = Min_dt_RKcontroller

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 51
nLevels = 2

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
