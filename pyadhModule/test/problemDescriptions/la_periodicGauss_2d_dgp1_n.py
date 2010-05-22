from pyadh import *
from pyadh.default_n import *
from la_periodicGauss_2d_p import *

timeOrder = 2
nStagesTime = timeOrder
#runCFL = 0.3
DT = None
nDTout = 1

timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller
usingSSPRKNewton = True

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#nn=41
nLevels = 1#4

subgridError = None

shockCapturing = None

if useHJ:
    numericalFluxType = HamiltonJacobi_DiagonalLesaintRaviart
else:
    numericalFluxType = Advection_DiagonalUpwind

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = SSPRKNewton #Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8
maxNonlinearIts =1

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
