from pyadh import *
from pyadh.default_n import *
from lad_shock_1d_p import *

#timeIntegration = BackwardEuler
#stepController = FixedStep
timeIntegration = FLCBDF
stepController = FLCBDF_controller
timeIntegration = BackwardEuler_cfl
#timeIntegration = ForwardEuler
#timeIntegration = ForwardEuler_A
stepController = MinDT_controller
runCFL = 0.099
#runCFL = None
nDTout = 10
DT=None
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)
#
#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn=101
nLevels = 1

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

massLumping=False

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False)

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#numericalFluxType = Advection_DiagonalUpwind
numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = numpy.array
multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

#conservativeFlux = {0:'pwl'}
