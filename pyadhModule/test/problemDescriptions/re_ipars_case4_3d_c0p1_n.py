from pyadh import *
from pyadh.default_n import *
from re_ipars_case4_3d_p import *

runCFL=None
timeIntegration = FLCBDF
stepController = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep

rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3

#timeIntegration = BackwardEuler
#stepController  = FixedStep
#systemStepControllerType = SplitOperator.Sequential_FixedStep
#DT = 1.0e-4#None#0.025#1.0e-1/timeScale
#nDTout =200#int(T/DT)

#tnList=[0.0,T]
#tnList=[0.0,1.0e-6,T]

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)

#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

triangleOptions = "pAYfen"
#triangleOptions = "pAfen"
#triangleOptions = "pqDAfen"

#nnx = 21
#nny = 21
#nnz = 21
nLevels = 1

massLumping=False

#subgridError = FFDarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=False)

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)

shockCapturing=None
shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True,nStepsToDelay=1)

numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts = 25

maxLineSearches = 25

fullNewtonFlag = True

tolFac = 1.0e-6

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = PETSc
multilevelLinearSolver = LU

levelLinearSolver = PETSc
levelLinearSolver = LU

linTolFac = 0.0001

#conservativeFlux = {0:'pwl',1:'pwl'}
