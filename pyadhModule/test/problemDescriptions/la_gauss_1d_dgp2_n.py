from pyadh import *
from pyadh.default_n import *
from la_gauss_1d_p import *

#mwf add
timeOrder =3
nStagesTime = timeOrder
runCFL = 0.1 #0.185
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
    

#mwf BackwardEuler 

timeIntegration = SSPRKwrap 
stepController=Min_dt_RKcontroller
DT = None
nDTout = 10

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)#SimplexLobattoQuadrature(nd,1)#

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 3
nLevels = 6

subgridError = None

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

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
