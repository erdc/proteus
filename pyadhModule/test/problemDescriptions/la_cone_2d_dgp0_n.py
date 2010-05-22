from pyadh import *
from pyadh.default_n import *
from la_cone_2d_p import *

timeOrder = 1
nStagesTime = timeOrder
runCFL = 0.45

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
#BackwardEuler, ForwardEuler
nStagesTime = 1
timeIntegration = BackwardEuler #SSPRKwrap
DT = None
nDTout = 10

femSpaces = {0:DG_Constants}

# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 3
nLevels = 4

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

#checkMass = True


archiveFlag = ArchiveFlags.EVERY_USER_STEP
