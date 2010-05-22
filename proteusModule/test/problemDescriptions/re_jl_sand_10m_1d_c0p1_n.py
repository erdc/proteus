from pyadh import *
from pyadh.default_n import *
from re_jl_sand_10m_1d_p import *

#timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
atol_u[0] = 1.0e-3
rtol_u[0] = 1.0e-3
timeIntegration = BackwardEuler
stepController  = FixedStep
systemStepControllerType = SplitOperator.Sequential_FixedStep
#DT = 1.0e-4#None#0.025#1.0e-1/timeScale
nDTout =100#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=101
#nLevels = 4
nn=101#3*2**2
nLevels = 1#10-2
#nn=3
#nLevels = 10
#tnList=[0.0,1.0e-8,T]
subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)

masslumping = False
#massLumping = True

shockCapturing = None
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.99,lag=True)
shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.99,lag=True,nStepsToDelay=1)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.99,lag=True)

numericalFluxType = None
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
#multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
multilevelNonlinearSolver = Newton

#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
#nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-4

maxNonlinearIts = 25#1001
maxLineSearches =25

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU
#multilevelLinearSolver = NI
#multilevelLinearSolver = Jacobi
#multilevelLinearSolver = GaussSeidel
#multilevelLinearSolver = StarILU

levelLinearSolver = LU
computeEigenvalues = False#True
#computeEigenvalues = False
#levelLinearSolver = LU#MGM

#linearSmoother = Jacobi
#linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.01

#conservativeFlux = {0:'pwl'}
