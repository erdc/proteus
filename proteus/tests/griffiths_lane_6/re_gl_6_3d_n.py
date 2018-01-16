from proteus import *
from proteus.default_n import *
from re_gl_6_3d_p import *

#unsteady
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-4
atol_u[0] = 1.0e-4
tnList = [0.0,1.0e-5,1.0e5]
timeIntegration = BackwardEuler
stepController = HeuristicNL_dt_controller
nonlinearIterationsFloor =5
nonlinearIterationsCeil=10
systemStepControllerType = SplitOperator.Sequential_MinModelStep
maxNonlinearIts=25
maxLineSearches=25
# maxNonlinearIts=100
# maxLineSearches=0#25
# #steady
# timeIntegration = PsiTCtte_new
# stepController = PsiTCtte_controller#SteadyStateIntegrator
# maxNonlinearIts=1
# maxLineSearches=0
#timeIntegration = NoIntegration
#stepController = Newton_controller
#maxNonlinearIts=100
#maxLineSearches=25
#systemStepControllerType = SplitOperator.Sequential_FixedStep
#tnList = [0.0,1.0]
# tnList = [0.0,1.0e-1]
#tnList=[float(i) for i in range(20)]

#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elementQuadrature = {('m',0):SimplexGaussQuadrature(nd,5),
                     'default':SimplexLobattoQuadrature(nd,1)}
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)


subgridError = None

massLumping = False

shockCapturing = None

numericalFluxType = None
#cek todo switch to sipg and adjust penalty
numericalFluxType = Richards_IIPG_exterior #need weak for parallel and global conservation 
#numericalFluxType = Richards_SIPG_exterior #need weak for parallel and global conservation 
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

multilevelNonlinearSolver = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.0

atol = 0.01*he#0.001*vFine

nl_rtol_res = 0.0
nl_atol_res = atol#0.001*vFine#1.0e-4
atol_res = {0:atol}
rtol_res = {0:0.0}
matrix = SparseMatrix

#multilevelLinearSolver =PETSc
#levelLinearSolver = PETSc
multilevelLinearSolver = LU
levelLinearSolver = LU
#multilevelLinearSolver =KSP_petsc4py
#levelLinearSolver = KSP_petsc4py

#conservativeFlux = {0:'pwl-bdm'}

conservativeFlux = None
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1
elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
archiveFlag = ArchiveFlags.EVERY_MODEL_STEP
