from pyadh import *
from pyadh.default_n import *
from sw_beach_erosion_board_2d_p import *

spaceOrder = 1; lagrangeBasis=True
if spaceOrder == 1:
    runCFL=0.1#0.3
elif spaceOrder == 2:
    runCFL=0.1
elif spaceOrder == 3:
    runCFL=0.05
else:
    runCFL=0.15

timeOrder = min(spaceOrder+1,3)
nStagesTime=timeOrder
timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller

nDTout = 101
tnList = [0.0 + i*T/(float(nDTout)-1.) for i in range(nDTout)]
archiveFlag = ArchiveFlags.EVERY_USER_STEP

if not lagrangeBasis and spaceOrder > 0:
    limiterType = DGlimiterPkMonomial2d
elif spaceOrder == 1:
    limiterType = DGlimiterDurlofskyP1Lagrange2d_Sw
    #limiterType = DGlimiterP1Lagrange2d
elif spaceOrder == 2:
    limiterType = DGlimiterP2Lagrange2d

stepController=Min_dt_RKcontroller
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

if spaceOrder == 1:
    if lagrangeBasis:
        femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis,
                     1:DG_AffineLinearOnSimplexWithNodalBasis,
                     2:DG_AffineLinearOnSimplexWithNodalBasis}
    else:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis,
                     1:DG_AffineP1_OnSimplexWithMonomialBasis,
                     2:DG_AffineP1_OnSimplexWithMonomialBasis}
elif spaceOrder == 2:
    if lagrangeBasis:
        femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis,
                     1:DG_AffineQuadraticOnSimplexWithNodalBasis,
                     2:DG_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis,
                     1:DG_AffineP2_OnSimplexWithMonomialBasis,
                     2:DG_AffineP2_OnSimplexWithMonomialBasis}
elif spaceOrder == 3:
    assert not lagrangeBasis
    femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis,
                 1:DG_AffineP3_OnSimplexWithMonomialBasis,
                 2:DG_AffineP3_OnSimplexWithMonomialBasis}
else:
    femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
                 1:DG_AffineP0_OnSimplexWithMonomialBasis,
                 2:DG_AffineP0_OnSimplexWithMonomialBasis}


if spaceOrder == 0:
    elementQuadrature = SimplexGaussQuadrature(nd,1)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)
else:
    elementQuadrature = SimplexGaussQuadrature(nd,2*spaceOrder+1)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)

nnx=51#2**(5-spaceOrder)+1
nny=21#2**(4-spaceOrder)+1

nLevels = 1

subgridError = ShallowWater_CFL(coefficients,nd,g)

massLumping=False

shockCapturing = None

numericalFluxType = RusanovNumericalFlux#ShallowWater_2D#RusanovNumericalFlux#

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
usingSSPRKNewton=True
#multilevelNonlinearSolver=Newton
levelNonlinearSolver = SSPRKNewton#Newton
#fullNewtonFlag = False

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001


