from proteus import *
from proteus.default_n import *
from vof_rotation_2d_p import *
from rotation2D import *

#multilevelNonlinearSolver  = NLNI
#levelNonlinearSolver = Newton
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = ExplicitConsistentMassMatrixForVOF

fullNewtonFlag = False
updateJacobian = False

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

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
    subgridError = None #Advection_ASGS(coefficients,nd,lag=False)
    shockCapturing = VOF3P.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_vof,lag=lag_shockCapturing_vof)
    if parallel or LevelModelType == VOF3P.LevelModel:
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
    elementQuadrature = CubeGaussQuadrature(nd,rotation_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,rotation_quad_order)
else:
    elementQuadrature = SimplexGaussQuadrature(nd,rotation_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,rotation_quad_order)

nonlinearSmoother = None#NLGaussSeidel

tolFac = 0.0
linTolFac = tolFac

nl_atol_res = 100*atolVolumeOfFluid
l_atol_res = atolVolumeOfFluid

maxNonlinearIts = 20
maxLineSearches = 0

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py#PETSc
    levelLinearSolver = KSP_petsc4py#PETSc
    linear_solver_options_prefix = 'vof_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

conservativeFlux = {}
if checkMass:
    auxiliaryVariables = [MassOverRegion()]
