from pyadh import *
from pyadh.default_n import *
from la_cone_2d_p import *

timeOrder = 3
nStagesTime = timeOrder
runCFL = 0.15
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

#BackwardEuler,SSPRKwrap
#timeIntegration = SSPRKwrap 
#stepController=Min_dt_RKcontroller


DT=None
nDTout = 20
timeOrder=1
nStagesTime=timeOrder
timeIntegration = BackwardEuler 
runCFL=0.5

#DT = 1.0e-3
#nDTout = T/DT

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 3
nLevels = 3

subgridError = None

massLumping = False

shockCapturing = None

numericalFluxType = Advection_DiagonalUpwind

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

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
