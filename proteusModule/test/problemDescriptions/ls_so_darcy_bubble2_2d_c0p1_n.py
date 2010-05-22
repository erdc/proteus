from pyadh import *
from pyadh.default_n import *
from ls_so_darcy_bubble2_2d_p import *

timeIntegration = BackwardEuler
#DT = 1.0e-6
#DT = T/10.0
DT = None #0.1
nDTout = int(T/0.1)
#timeIntegration = ForwardEuler_H
#runCFL = 0.1
#DT = 1.0e-5
#timeIntegrator = ForwardIntegrator


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 5
nDTout = 100

#subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd)

massLumping = False

numericalFluxType = None

#shockCapturing = None
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
