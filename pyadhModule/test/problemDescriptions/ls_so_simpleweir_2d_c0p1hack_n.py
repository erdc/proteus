from pyadh import *
from pyadh.default_n import *
from ls_so_simpleweir_2d_p import *

timeIntegration = BackwardEuler
timeIntegration = ForwardEuler_H
timeIntegrator = ForwardIntegrator

timeIntegrator = testStuff.ExplicitLevelSetIntegrator

runCFL = 0.1

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=3
nLevels = 5
DT = None
nDTout = 100

subgridError = HamiltonJacobi_ASGS(coefficients,nd)

massLumping = False

numericalFluxType = None

shockCapturing = ResGrad_SC(coefficients,nd)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

maxNonlinearIts = 1000

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
