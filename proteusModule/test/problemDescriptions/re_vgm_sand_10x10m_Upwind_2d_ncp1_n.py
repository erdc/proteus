from pyadh import *
from pyadh.default_n import *
from re_vgm_sand_10x10m_Upwind_2d_p import *

#timeIntegration = BackwardEuler
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3
DT = None
femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nnx=51#3
nny=51
#nn = 51
#mwf something wrong with 61?
nLevels = 1#1
nDTout = 1#int(T/DT)

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-7

maxNonlinearIts = 20#0
maxLineSearches = 5

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'p1-nc'}
