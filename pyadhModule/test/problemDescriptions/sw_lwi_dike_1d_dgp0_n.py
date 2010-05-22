from pyadh import *
from pyadh.default_n import *
from sw_lwi_dike_1d_p import *

runCFL=0.33
#runCFL=0.15
timeOrder = 1
timeIntegration = LinearSSPRKintegration
stepController=Min_dt_RKcontroller
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

#nDTout=1001
femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
             1:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
#elementQuadrature = SimplexGaussQuadrature(nd,2)
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

nn=201
nLevels = 1

subgridError = ShallowWater_CFL(coefficients,nd,g)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

massLumping=False

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#numericalFluxType = Advection_DiagonalUpwind
numericalFluxType = ShallowWater_1D

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
