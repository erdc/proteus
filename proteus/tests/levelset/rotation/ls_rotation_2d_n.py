from proteus import *
from proteus.default_n import *
try:
    from .ls_rotation_2d_p import *
    from .rotation2D import *
except:
    from ls_rotation_2d_p import *
    from rotation2D import *
nd = 2

if timeIntegration_ls == "be":
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
elif timeIntegration_ls == "vbdf":
    timeIntegration = VBDF
    stepController = Min_dt_cfl_controller
elif timeIntegration_ls == "flcbdf":
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
elif timeIntegration_ls == "rk":
    if cDegree_ls == -1:
        timeIntegration = LinearSSPRKPIintegration
    else:
        timeIntegration = LinearSSPRKintegration        
    stepController=Min_dt_RKcontroller
    timeOrder = pDegree_ls+1
    nStagesTime = timeOrder
else:
    raise RuntimeError

if useHex:
    if pDegree_ls == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_ls == 2:
        femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}#this is hardwired to p2 right now
    else:
        print("pDegree_ls = %s not recognized " % pDegree_ls)
    elementQuadrature = CubeGaussQuadrature(nd,rotation_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,rotation_quad_order)
else:
    if pDegree_ls == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ls == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print("pDegree_ls = %s not recognized " % pDegree_ls)
    elementQuadrature = SimplexGaussQuadrature(nd,rotation_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,rotation_quad_order)

subgridError = None

subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=True)

massLumping = False

numericalFluxType = None

shockCapturing = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_ls,lag=lag_shockCapturing_ls)

numericalFluxType = DoNothing


multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton

nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
nonlinearSmoother = None

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolLevelSet
l_atol_res = atolLevelSet

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py#PETSc
    levelLinearSolver = KSP_petsc4py#PETSc
    linear_solver_options_prefix = 'ncls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU

conservativeFlux = {}
maxNonlinearIts = 5
maxLineSearches = 0

#checkMass = True

if not applyCorrection and checkMass:
   auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryLS("rotation2dnc"+repr(lRefinement))]
