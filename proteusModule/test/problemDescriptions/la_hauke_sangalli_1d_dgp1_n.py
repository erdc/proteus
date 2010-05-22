from pyadh import *
from pyadh.default_n import *
from la_hauke_sangalli_1d_p import *

RKtimeIntegration = False
if RKtimeIntegration:
    timeOrder = 2
    nStagesTime = timeOrder
    timeIntegration = LinearSSPRKintegration
    stepController  = Min_dt_RKcontroller
    runCFL = 0.1
else:
    stepController = Min_dt_controller
    nn = haukeSangalli_nn
    timeIntegration = BackwardEuler
    DT = haukeSangalli_dt
    nDTout = int(T/DT)
 
nLevels = 1
nn = haukeSangalli_nn
nDTout = int(T/haukeSangalli_dt)

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)


subgridError = None
if RKtimeIntegration:
    numericalFluxType = Advection_DiagonalUpwind#_Diffusion_IIPG
else:
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
multilevelNonlinearSolver  = NLNI

maxNonlinearIts =1

maxLineSearches =0

levelNonlinearSolver = Newton

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
