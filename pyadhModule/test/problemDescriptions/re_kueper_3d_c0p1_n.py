from pyadh import *
from pyadh.default_n import *
from re_kueper_3d_p import *

# timeIntegration = BackwardEuler
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-2
atol_u[0] = 1.0e-2
tnList=[0.0,T]

#DT = 1.0e-6/timeScale
#1.0e-3/timeScale
#nDTout = 100
# timeIntegration = BackwardEuler
# stepController  = FixedStep
# systemStepControllerType = SplitOperator.Sequential_FixedStep
# #DT = 1.0e-4#None#0.025#1.0e-1/timeScale
# nDTout =10#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

#appears to work
elementQuadrature = {'default':SimplexLobattoQuadrature(nd,1),
                     ('u',0):SimplexGaussQuadrature(nd,3)}
elementBoundaryQuadrature = {'default':SimplexLobattoQuadrature(nd-1,1),
                             ('u',0):SimplexGaussQuadrature(nd-1,3)}

elementQuadrature = SimplexLobattoQuadrature(nd,1)

elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nLevels = 1
#triangleOptions="pAq30Dena%f" % ((0.5*top*right)*0.0005,)
#triangleOptions="pAq30Dena%f" % ((0.5*top*right)*0.001,)
triangleOptions = "pAfena%e" % ((top*right*top*0.1/6.0)*0.001,)
subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

shockCapturing = None
shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.99,lag=True,nStepsToDelay=1)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)

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

maxNonlinearIts = 25
maxLineSearches =25

matrix = SparseMatrix
#matrix = numpy.array

multilevelLinearSolver = LU
multilevelLinearSolver = PETSc
#multilevelLinearSolver = NI

levelLinearSolver = LU
levelLinearSolver = PETSc
#levelLinearSolver = MGM

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'pwl'}
conservativeFlux = {0:'point-eval'}
