from proteus import *
from proteus.default_n import *
from .vof_vortex_3d_p import *
from .vortex import *

multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton
fullNewtonFlag = True

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
    limiterType =   {0:TimeIntegration.DGlimiterPkMonomial2d}
else:
    raise RuntimeError

if cDegree_vof==0:
    if useHex:
        if pDegree_vof==1:
            femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
        elif pDegree_vof==2:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        if pDegree_vof==1:
            femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
        elif pDegree_vof==2:
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    subgridError = Advection_ASGS(coefficients,nd,lag=False)
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=shockCapturingFactor_vof,lag=True)
    if parallel or LevelModelType == VOF.LevelModel:
        numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

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
    limiterType =   {0:TimeIntegration.DGlimiterPkMonomial2d}

if useHex:
    elementQuadrature = CubeGaussQuadrature(nd,vortex_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,vortex_quad_order)
else:
    elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nonlinearSmoother = NLGaussSeidel

tolFac = 0.0

nl_atol_res = atolVolumeOfFluid

maxNonlinearIts = 50

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    linear_solver_options_prefix = 'vof_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = {}
