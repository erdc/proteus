from pyadh import *
from pyadh.default_n import *
from ls_vortex_2d_p import *
from vortex import *

timeOrder = 2
nStagesTime = timeOrder
DT = None
nDTout = 1

timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller
usingSSPRKNewton = True


femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

if useHJ:
    numericalFluxType = HamiltonJacobi_DiagonalLesaintRaviart
else:
    numericalFluxType = Advection_DiagonalUpwind


subgridError = None

massLumping = False

shockCapturing = None


multilevelNonlinearSolver  = Newton

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
