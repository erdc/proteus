from pyadh import *
from pyadh.default_n import *
from lwi_dike import *
from ls_lwi_dike_2d_p import *


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
#timeIntegration = ForwardEuler
timeIntegration = BackwardEuler
#runCFL = 1.5
#timeIntegration = FLCBDF
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



femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,lwi_dike_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,lwi_dike_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind

shockCapturing = None

#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 0.01

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
