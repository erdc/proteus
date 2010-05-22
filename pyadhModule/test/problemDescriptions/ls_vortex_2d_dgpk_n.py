from pyadh import *
from pyadh.default_n import *
from ls_vortex_2d_p import *


timeOrder = 4
nStagesTime = timeOrder
timeIntegration = LinearSSPRKintegration
stepController=Min_dt_RKcontroller
usingSSPRKNewton = True

DT = None
nDTout = 1
femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,6)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


subgridError = None

massLumping  = False

shockCapturing = None

if useHJ:
    numericalFluxType = HamiltonJacobi_DiagonalLesaintRaviart
else:
    numericalFluxType = Advection_DiagonalUpwind


multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = SSPRKNewton#Newton

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
