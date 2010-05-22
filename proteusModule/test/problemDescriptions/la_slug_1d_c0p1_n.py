from pyadh import *
from pyadh.default_n import *
from la_slug_1d_p import *

nLevels = 1
#BackwardEuler
stepController = Min_dt_controller
if harariTestProblem:
    timeIntegration = BackwardEuler
    DT = L[0]/(nn-1.)*L[0]/(nn-1.)*0.0025/a0
    nDTout = int(T/DT)
    nn = 6
else:
    timeIntegration = BackwardEuler#BackwardEuler_cfl
    nn = 51
    runCFL = 0.1#0.1
    DT = runCFL/b0/(nn-1.)
    nDTout = int(T/DT)

#timeIntegration = FLCBDF
#stepController = FLCBDF_controller
#atol_u[0] = 1.0e-2
#rtol_u[0] = 1.0e-2
#timeOrder= 5

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)


subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='1',lag=False)
#subgridError = AdvectionDiffusionReactionHaukeSangalliInterpolant_ASGS(coefficients,nd,stabFlag='1',lag=True,interpolationFemSpaceType = femSpaces[0])
#subgridError = AdvectionDiffusionReactionHaukeSangalliInterpolantWithTransientSubScales_ASGS(coefficients,nd,stabFlag='1',lag=True,interpolationFemSpaceType = femSpaces[0],
#                                                                                             trackSubScales=True)
subgridError = AdvectionDiffusionReactionTransientSubscales_ASGS(coefficients,nd,stabFlag='1',lag=False,trackSubScales=True,useHarariDirectly=harariTestProblem)
numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1)

multilevelNonlinearSolver  = NLNI

maxNonlinearIts =10#5

maxLineSearches =0

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

