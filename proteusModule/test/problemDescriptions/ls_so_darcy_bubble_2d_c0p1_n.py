from pyadh import *
from pyadh.default_n import *
from ls_so_darcy_bubble_2d_p import *

timeIntegration = BackwardEuler
#DT = 1.0e-6
#DT = T/10.0
DT = 0.01
nDTout = int(T/DT)
#timeIntegration = ForwardEuler_H
#runCFL = 0.1
#DT = 1.0e-5
#timeIntegrator = ForwardIntegrator


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
nLevels = 5


#subgridError = None
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)
subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True) #doesn't work now for some reason

massLumping = False

numericalFluxType = None

#shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd,lag=True)

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
