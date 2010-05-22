from pyadh import *
from pyadh.default_n import *
from re_forsyth2_ss_2d_p import *

timeIntegration = NoIntegration

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = {'default':SimplexGaussQuadrature(nd,1),
#                     ('u',0):SimplexGaussQuadrature(nd,1)}
#elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,1),
#                             ('u',0):SimplexGaussQuadrature(nd-1,1)}

#should be ok with either post processing, since K is diagonal, no mass terms
elementQuadrature = {'default':SimplexGaussQuadrature(nd,2),
                     ('u',0):SimplexGaussQuadrature(nd,2)}
#can insert nodal quadrature here, or let SimTools compute. Inserting here can change convergence
#at least for stabilized method I believe
#                     'NodalQuadrature':SimplexLobattoQuadrature(nd,1)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,2),
                             ('u',0):SimplexGaussQuadrature(nd-1,2)}

nn=3
nLevels = 3
triangleOptions+="A"# for element type data

subgridError = None

massLumping = False

numericalFluxType = None

shockCapturing = None

multilevelNonlinearSolver = NLNI #Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 1.0e-6

nl_atol_res = 1.0e-6

maxNonlinearIts = 100#1

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'p1-nc'}
