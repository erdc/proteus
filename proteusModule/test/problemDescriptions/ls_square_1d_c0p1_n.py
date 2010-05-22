from pyadh import *
from pyadh.default_n import *
from ls_square_1d_p import *
from square import *

if timeIntegration_ls == "BE":
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
elif timeIntegration_ls == "FLCBDF":
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller
elif timeIntegration_ls == "RK":
    if cDegree_ls == -1:
        timeIntegration = LinearSSPRKPIintegration
    else:
        timeIntegration = LinearSSPRKintegration        
    stepController=Min_dt_RKcontroller
    timeOrder = pDegree_ls+1
    nStagesTime = timeOrder
else:
    raise RuntimeError

if cDegree_ls==0:
    if pDegree_ls == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ls ==2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    if LevelModelType == NCLS.OneLevelNCLS:
        subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=True)
    else:
        subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=shockCapturingFactor_ls,lag=True)
    if LevelModelType == NCLS.OneLevelNCLS:
        numericalFluxType = DoNothing
        
elif cDegree_ls==-1:
    if pDegree_ls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    elif pDegree_ls==3:
        femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
    numericalFluxType = HamiltonJacobi_DiagonalLesaintRaviart
    #limiterType =   {0:TimeIntegration.DGlimiterPkMonomial1d}
    limiterType =   {0:None}
    subGridError = None

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)



massLumping = False


#shockCapturing = ResGrad_SC(coefficients,nd)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False)
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

