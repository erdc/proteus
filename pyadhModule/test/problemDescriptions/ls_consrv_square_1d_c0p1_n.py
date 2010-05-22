from pyadh import *
from pyadh.default_n import *
from ls_consrv_square_1d_p import *
from square import *


#timeIntegrator = ForwardIntegrator
#timeIntegration = NoIntegration
timeIntegration = NoIntegration
stepController = Newton_controller

if cDegree_ls==0:
    if pDegree_ls==1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ls==2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
elif cDegree_ls==-1:
    if pDegree_ls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    numericalFluxType = Diffusion_SIPG

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False


shockCapturing = None


#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolConservation

maxNonlinearIts = 100

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 1.0e-6

conservativeFlux = None
