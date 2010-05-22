from pyadh import *
from pyadh.default_n import *
from ls_square_1d_p import *
from square import *

timeOrder =3
nStagesTime = timeOrder

DT=None
#now in suqare
#runCFL = 0.185

class SSPRKwrap(SSPRKintegration):
    """
    wrap SSPRK so default constructor uses
    the order and runCFL without
    changing Transport
    """
    def __init__(self,vt):
        SSPRKintegration.__init__(self,vt,timeOrder,runCFL,usingSSPRKNewton=False)
        return
    #
#end wrapper
    
timeIntegration = SSPRKwrap 
stepController=Min_dt_RKcontroller
nDTout = 1

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

#now in square
#nn = 17#
#nLevels = 1#6

subgridError = None

massLumping = False

numericalFluxType = Advection_DiagonalUpwind

shockCapturing = None

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

#checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
