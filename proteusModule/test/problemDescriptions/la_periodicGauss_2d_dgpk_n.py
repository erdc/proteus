from pyadh import *
from pyadh.default_n import *
from la_periodicGauss_2d_p import *

timeOrder = 4
nStagesTime = timeOrder
runCFL = 0.125
DT = None
nDTout = 1

timeIntegration = LinearSSPRKintegration
stepController=Min_dt_RKcontroller

femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,6)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,6)

#nn=11
nLevels = 1#4

subgridError = None

shockCapturing = None

if useHJ:
    numericalFluxType = HamiltonJacobi_DiagonalLesaintRaviart
else:
    numericalFluxType = Advection_DiagonalUpwind
    
multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

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
