from pyadh import *
from pyadh.default_n import *
from re_vgm_sand_10x10m_l2_ss_2d_p import *

timeIntegration = NoIntegration

DT = 1.0

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}


#if use L2 proj still need higher order quad to get Jacobian right,
#still should get back same integrals since stiffness shape/test
#function terms are all integrated exactly with barycenter quadrature
#(they're constant)
elementQuadrature = {'default':SimplexGaussQuadrature(nd,2),
                     ('u',0):SimplexGaussQuadrature(nd,2)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,2),
                             ('u',0):SimplexGaussQuadrature(nd-1,2)}

nn=3
nLevels = 4

subgridError = None

massLumping = False

numericalFluxType = None

shockCapturing = None

multilevelNonlinearSolver = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

maxNonlinearIts = 100#1

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'p1-nc'}
