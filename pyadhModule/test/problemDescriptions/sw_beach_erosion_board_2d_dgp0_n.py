from pyadh import *
from pyadh.default_n import *
from sw_beach_erosion_board_2d_p import *

#runCFL=0.33
runCFL=0.15
timeOrder = 1
nStagesTime = 1
timeIntegration = LinearSSPRKintegration
stepController=Min_dt_RKcontroller
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

#nDTout=1001
femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
             1:DG_AffineP0_OnSimplexWithMonomialBasis,
             2:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
#elementQuadrature = SimplexGaussQuadrature(nd,2)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

nnx=101
nny=21
nLevels = 1

subgridError = ShallowWater_CFL(coefficients,nd,g)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

massLumping=False

shockCapturing = None

numericalFluxType = ShallowWater_2D

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
usingSSPRKNewton=True
multilevelNonlinearSolver=Newton
levelNonlinearSolver = SSPRKNewton#Newton
#fullNewtonFlag = False

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

#conservativeFlux = {0:'pwl'}

archiveFlag = ArchiveFlags.EVERY_USER_STEP
