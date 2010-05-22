from pyadh import *
from pyadh.default_n import *
from twp_stokes_ls_fc_outstep_2d_p import *

timeIntegration = BackwardEuler

DT=T/100.0
nDTout=100

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineQuadraticOnSimplexWithNodalBasis,
             3:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 3
nLevels = 4

subgridError = TwophaseStokes_LS_FC_ASGS(coefficients,nd)

shockCapturing = None#ResGrad_SC(coefficients,nd)

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
