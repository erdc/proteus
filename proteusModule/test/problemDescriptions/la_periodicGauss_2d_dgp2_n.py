from pyadh import *
from pyadh.default_n import *
from la_periodicGauss_2d_p import *

timeOrder = 3
nStagesTime = timeOrder
#runCFL = 0.185
DT = None
nDTout = 1

#BackwardEuler,SSPRKwrap
timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#nn=21
nLevels = 1#4

subgridError = None

shockCapturing = None

if useHJ:
    numericalFluxType = HamiltonJacobi_DiagonalLesaintRaviart
else:
    numericalFluxType = Advection_DiagonalUpwind
    
multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = SSPRKNewton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.00

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
