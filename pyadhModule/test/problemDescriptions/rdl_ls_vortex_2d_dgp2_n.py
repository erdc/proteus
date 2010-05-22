from pyadh import *
from pyadh.default_n import *
from rdl_ls_vortex_2d_p import *
from vortex import *

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
timeIntegration = BackwardEuler
#timeIntegration = OuterTheta
#timeIntegration = FLCBDF
timeIntegration = ForwardEuler_A
timeIntegrator = ForwardIntegrator
timeIntegration = SSPRKwrap 
stepController=Min_dt_RKcontroller
#SSPRK
timeOrder = 2
nStagesTime = timeOrder

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)


subgridError = None
#subgridError = Advection_ASGS(coefficients,nd,lag=False)
#subgridError = Advection_ASGS(coefficients,nd,lag=True)

massLumping = False

#numericalFluxType = None
numericalFluxType = Advection_DiagonalUpwind

shockCapturing = None#ResGrad_SC(coefficients,nd)
#shockCapturing = ResGrad_SC(coefficients,nd)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)
multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 0.01/(nn -1.0 )

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None#{0:'dg-point-eval'}

#checkMass = True


archiveFlag = ArchiveFlags.EVERY_USER_STEP
