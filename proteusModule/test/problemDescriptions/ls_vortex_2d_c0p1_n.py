from pyadh import *
from pyadh.default_n import *
from ls_vortex_2d_p import *
from vortex import *


if timeIntegration_ls == "BE":
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    if timeOrder == 2:
        timeIntegration = VBDF
        stepController  = Min_dt_cfl_controller
    
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

if pDegree_ls == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elif pDegree_ls == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
else:
    print "pDegree_ls = %s not recognized " % pDegree_ls

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)


subgridError = None
if useHJ:
    if LevelModelType == NCLS.OneLevelNCLS:
        subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=True)
    else:
        subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)
else:
    subgridError = Advection_ASGS(coefficients,nd,lag=True)

massLumping = False

numericalFluxType = None

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=shockCapturingFactor_ls,lag=True)
if LevelModelType == NCLS.OneLevelNCLS or parallel:
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
    limiterType =   {0:TimeIntegration.DGlimiterPkMonomial2d}
    subGridError = None
    assert not parallel, "DG requires serial for now"

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolLevelSet

maxNonlinearIts = 2

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {}


if not applyCorrection:
   auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryLS("vortex2dnc"+`lRefinement`)]
