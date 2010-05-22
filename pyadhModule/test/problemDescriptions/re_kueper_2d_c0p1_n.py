from pyadh import *
from pyadh.default_n import *
from re_kueper_2d_p import *

#timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3
tnList=[0.0,T]
#DT = 1.0e-6/timeScale
#1.0e-3/timeScale
#nDTout = 100
# timeIntegration = BackwardEuler
# stepController  = FixedStep
# systemStepControllerType = SplitOperator.Sequential_FixedStep
# #DT = 1.0e-4#None#0.025#1.0e-1/timeScale
# nDTout =100#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

#appears to work
#elementQuadrature = SimplexLobattoQuadrature(nd,1)

#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nLevels = 1
#triangleOptions="pAq30Dena%f" % ((0.5*top*right)*0.0005,)
triangleOptions="pAq30Dena%f" % ((0.5*top*right)*0.001,)

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

shockCapturing = None
shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True,nStepsToDelay=1)

numericalFluxType = None
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

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

maxNonlinearIts = 10
maxLineSearches =10

matrix = SparseMatrix
#matrix = numpy.array

multilevelLinearSolver = LU
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = MGM

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'pwl'}

