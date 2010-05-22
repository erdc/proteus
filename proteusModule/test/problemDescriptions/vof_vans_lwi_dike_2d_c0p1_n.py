from pyadh import *
from pyadh.default_n import *
from lwi_dike import *
from vof_vans_lwi_dike_2d_p import *


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
timeIntegrator = ForwardIntegrator
timeIntegration = BackwardEuler
timeOrder = 1
nStagesTime = timeOrder
#timeIntegration = BackwardEuler
#BackwardEuler,SSPRKwrap
#timeIntegration = SSPRKwrap 
#stepController=Min_dt_RKcontroller
#SSPRK
#timeOrder = 3
#nStagesTime = timeOrder
#runCFL = 0.15
#DT=1.e-5
#DT=None
#nDTout = 100

femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,lwi_dike_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,lwi_dike_quad_order)

subgridError = None

massLumping = False

numericalFluxType = Advection_DiagonalUpwind

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 1.0e-4

nl_atol_res = 1.0e-4

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
