from pyadh import *
from pyadh.default_n import *
from sw_hump_1d_p import *

timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[1] = 1.0e-4
rtol_u[2] = 1.0e-4
atol_u[1] = 1.0e-4
atol_u[2] = 1.0e-4
#timeIntegration=BackwardEuler
timeIntegration=ForwardEuler
stepController=FixedStep

runCFL=0.45
timeOrder = 3
class SSPRKwrap(LinearSSPRKintegration):
    def __init__(self,vt):
        LinearSSPRKintegration.__init__(self,vt,timeOrder,runCFL)
        return
timeIntegration = SSPRKwrap 
stepController=Min_dt_RKcontroller
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

nDTout=1001
femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
             1:DG_AffineP0_OnSimplexWithMonomialBasis}
# femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis,
#             1:DG_AffineP1_OnSimplexWithMonomialBasis}
#femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#             1:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#             1:C0_AffineQuadraticOnSimplexWithNodalBasis}

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn=101
nLevels = 1

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

massLumping=False

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#numericalFluxType = Advection_DiagonalUpwind
numericalFluxType = ShallowWater_1D

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
#levelNonlinearSolver = testStuff.SSPRKNewton#Newton

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = numpy.array
multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

#conservativeFlux = {0:'pwl'}
