from pyadh import *
from pyadh.default_n import *
from re_bcb_sand_10x10m_2d_p import *

#timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
atol_u[0] = 1.0e-4
rtol_u[0] = 1.0e-4
tnList = [0.0,1.0e-8,T]
#DT = None#0.025#1.0e-1/timeScale
#nDTout = 100#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,4)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nnx=61
nny=61
nLevels = 1

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)

massLumping = False

numericalFluxType = None
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

shockCapturing = None
#shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.99,lag=True,nStepsToDelay=10)

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
#multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
multilevelNonlinearSolver = Newton

#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

#nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

maxNonlinearIts = 10#1001
maxLineSearches =10

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = PETSc
#levelLinearSolver = MGM

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

#conservativeFlux = {0:'pwl-bdm'}
parallelPartitioningType = MeshParallelPartitioningTypes.element
#default number of layers to use > 1 with element partition means
#C0P1 methods don't need to do communication in global element assembly
#nodal partitioning does not need communication for C0P1 (has overlap 1) regardless
nLayersOfOverlapForParallel = 1
