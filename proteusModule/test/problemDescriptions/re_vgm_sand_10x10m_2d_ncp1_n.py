from pyadh import *
from pyadh.default_n import *
from re_vgm_sand_10x10m_2d_p import *

timeIntegration = BackwardEuler

runCFL = 0.9

DT = 1.0e-3/timeScale
nDTout = 100
femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

#appears to work, make sure using P1 post processing that
#allows for "full" mass integral
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

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
#multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
multilevelNonlinearSolver = Newton

#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

#nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

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

conservativeFlux = {0:'p1-nc'}
