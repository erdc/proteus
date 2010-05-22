from pyadh import *
from pyadh.default_n import *
from re_forsyth2_ss_2d_p import *

timeIntegration = NoIntegration
#timeIntegration = PsiTCtte
#timeIntegrator = SteadyStateIntegrator

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#appears to work
#elementQuadrature = {'default':SimplexLobattoQuadrature(nd,3),
#                     ('u',0):SimplexGaussQuadrature(nd,3)}
#have to have higher order quadrature on element boundary to get velPPv2 with
#neumann bc's enforced explicitly
#elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,2),
#                             ('u',0):SimplexGaussQuadrature(nd-1,2)}

#mwf looks like works too only with stabilization?
elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)


nn=3
nLevels = 5

#subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

numericalFluxType = None

#shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
#multilevelNonlinearSolver = Newton

#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

#nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 1.0e-6 #1.0e-6

nl_atol_res = 1.0e-6 #1.0e-6

maxNonlinearIts = 100#1

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = MGM

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'pwl'}#{0:'point-eval'}#{0:'sun-gs-rt0'}#{0:'sun-rt0'}#{0:'pwl'}

#if using sun-gs-rt0 need weak dirichlet conditions
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

