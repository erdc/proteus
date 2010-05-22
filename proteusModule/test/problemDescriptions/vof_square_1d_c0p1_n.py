from pyadh import *
from pyadh.default_n import *
from vof_square_1d_p import *
from square import *

if timeIntegration_vof == "BE":
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
elif timeIntegration_vof == "FLCBDF":
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller
elif timeIntegration_vof == "RK":
    if cDegree_vof == -1:
        timeIntegration = LinearSSPRKPIintegration
    else:
        timeIntegration = LinearSSPRKintegration        
    stepController=Min_dt_RKcontroller
    #timeIntegration = ForwardEuler
    #stepController=Min_dt_controller
    timeOrder = pDegree_vof+1
    nStagesTime = timeOrder
    limiterType =   {0:TimeIntegration.DGlimiterPkMonomial1d}
else:
    raise RuntimeError

if cDegree_vof==0:
    if pDegree_vof==1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_vof==2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    subgridError = Advection_ASGS(coefficients,nd,lag=True)
    #shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=shockCapturingFactor_vof,lag=True)
elif cDegree_vof==-1:
    if pDegree_vof==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_vof==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_vof==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    elif pDegree_vof==3:
        femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
    numericalFluxType = Advection_DiagonalUpwind
    subGridError = None

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

# elementQuadrature = SimplexLobattoQuadrature(nd,vortex_quad_order)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,vortex_quad_order)

massLumping = False

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {}

