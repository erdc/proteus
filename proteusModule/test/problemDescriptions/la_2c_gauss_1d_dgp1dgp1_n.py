from pyadh import *
from pyadh.default_n import *
from la_2c_gauss_1d_p import *

timeOrder = 3
nStagesTime = timeOrder
runCFL=0.1
class SSPRKwrap(SSPRKintegration):
    """
    wrap SSPRK so default constructor uses
    the order I want and runCFL without
    changing VectorTransport
    """
    def __init__(self,vt):
        SSPRKintegration.__init__(self,vt,timeOrder,runCFL)
        return
    #
#end wrapper

timeIntegration = SSPRKwrap 
stepController=Min_dt_RKcontroller

DT = None
nDTout = 10


femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis,
             1:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 3
nLevels = 5

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




archiveFlag = ArchiveFlags.EVERY_USER_STEP
