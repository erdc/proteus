from proteus import *
from proteus.default_n import *
from re_vgm_sand_10x10x10_3d_p import *

#timeIntegration = BackwardEuler
#timeIntegration = FLCBDF
#stepController = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3
timeIntegration = BackwardEuler
stepController = HeuristicNL_dt_controller#FixedStep
#nDTout = 1#int(T/DT)#int(T/DT) #100#int(T/DT)
#for controlling time stepping
nonlinearIterationsFloor = 4
nonlinearIterationsCeil  = 8
dtNLgrowFactor = 2
dtNLreduceFactor = 0.5
dtNLfailureReduceFactor = 0.5
useInitialGuessPredictor= True
stepExact = True
nDTout = 100
DT = T/nDTout 
tnList = [0.0,min(1.0e-5,DT/2.0)]+[i*DT for i  in range(1,nDTout+1)]

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
#pwl-bdm requires more accurate quadrature on element boundaries
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
#pwl-bdm requires more accurate quadrature on element boundaries
#nnx=61
#nny=4
#nnz=61
nnx=21
nny=21
nnz=21
nLevels =1

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)

massLumping = False

numericalFluxType = None
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

shockCapturing = None
#shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.75,lag=True,nStepsToDelay=1)

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

nl_atol_res = 1.0e-6

maxNonlinearIts = 10
maxLineSearches = 25
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
parallelPartitioningType = MeshParallelPartitioningTypes.node
#default number of layers to use > 1 with element partition means
#C0P1 methods don't need to do communication in global element assembly
#nodal partitioning does not need communication for C0P1 (has overlap 1) regardless
nLayersOfOverlapForParallel = 1
