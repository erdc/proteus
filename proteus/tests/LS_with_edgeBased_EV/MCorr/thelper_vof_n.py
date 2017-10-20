from proteus import *
from proteus.default_n import *
from thelper_vof_p import *
from thelper_cons_ls import *

multilevelNonlinearSolver  = Newton
if ct.STABILIZATION_TYPE_vof == 0: #SUPG
    levelNonlinearSolver = Newton
    fullNewtonFlag = True
    updateJacobian = True
    timeIntegration = BackwardEuler_cfl
else:
    fullNewtonFlag = False
    updateJacobian = False
    timeIntegration = VOF.RKEV # SSP
    levelNonlinearSolver = ExplicitConsistentMassMatrixForVOF

stepController = Min_dt_controller
timeOrder = ct.SSPOrder
nStagesTime = ct.SSPOrder

if useHex:
    hex=True
    quad=True
    if pDegree_vof==1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_vof==2:
        if useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_vof
    elementQuadrature = CubeGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,quad_order)
else:
    if pDegree_vof == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_vof == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_vof
    elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)
    
if parallel or LevelModelType == VOF.LevelModel:
    numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior
    
#subgridError = Advection_ASGS(coefficients,nd,lag=False)
shockCapturing = VOF.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_vof,lag=lag_shockCapturing_vof)

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'vof_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU


