from pyadh import *
from pyadh.default_n import *
from sw_hump_2d_p import *

implicit=False

if implicit:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller
    rtol_u[0] = 1.0e-4
    rtol_u[1] = 1.0e-4
    rtol_u[2] = 1.0e-4
    atol_u[0] = 1.0e-4
    atol_u[1] = 1.0e-4
    atol_u[2] = 1.0e-4
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis}
    elementQuadrature = SimplexLobattoQuadrature(nd,1)
    elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
    multilevelNonlinearSolver  = Newton
    
    levelNonlinearSolver = Newton

    fullNewtonFlag = True
else:
    runCFL=0.25
    timeOrder = 1
    class SSPRKwrap(LinearSSPRKintegration):
        def __init__(self,vt):
            LinearSSPRKintegration.__init__(self,vt,timeOrder,runCFL)
            return
    timeIntegration = SSPRKwrap 
    stepController=Min_dt_RKcontroller
    systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
    #nDTout=101
    femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
                 1:DG_AffineP0_OnSimplexWithMonomialBasis,
                 2:DG_AffineP0_OnSimplexWithMonomialBasis}
#    femSpaces = {0:DG_Constants,
#                 1:DG_Constants,
#                 2:DG_Constants}
    elementQuadrature = SimplexGaussQuadrature(nd,1)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)
    numericalFluxType = ShallowWater_2D
    multilevelNonlinearSolver  = Newton
    levelNonlinearSolver = testStuff.SSPRKNewton
    fullNewtonFlag = True
#    femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis,
#                 1:DG_AffineP1_OnSimplexWithMonomialBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#             1:C0_AffineQuadraticOnSimplexWithNodalBasis}

nnx=81
nny=81
nLevels = 1

subgridError = ShallowWater_CFL(coefficients,nd,g)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

massLumping=False

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#numericalFluxType = Advection_DiagonalUpwind


tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = numpy.array
multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

#conservativeFlux = {0:'pwl'}
