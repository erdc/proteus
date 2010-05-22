from pyadh import *
from pyadh.default_n import *
from redist_circle_1d_p import *


#timeIntegration = NoIntegration
#stepController = Newton_controller
timeIntegration = BackwardEuler
stepController = Osher_controller
#timeIntegration = PsiTCtte
#stepController = PsiTCtte_controller

runCFL = 0.1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}


elementQuadrature = SimplexGaussQuadrature(nd,3)
#
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=21
nLevels = 1#5
#mwf force a small initial dt?
DT = 1.0
nDTout = 1

#subgridError = HamiltonJacobi_ASGS(coefficients,nd)
subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)

massLumping = False

#shockCapturing = ResGrad_SC(coefficients,nd)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.0,lag=True)
numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#these need to be set appropriately for pseudo-transient
tolFac = 0.1
nl_atol_res = 1.0e-3

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
