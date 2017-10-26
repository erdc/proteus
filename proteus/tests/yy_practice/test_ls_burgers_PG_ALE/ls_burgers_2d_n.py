from proteus import *
from proteus.default_n import *
from ls_burgers_2d_p import *
from burgers2D import *

nd = 2

quad = False
if ct.useHex:
    quad = True
    if pDegree_ls == 1:
        femSpaces = {0: C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_ls == 2:
        # this is hardwired to p2 right now
        femSpaces = {0: C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print "pDegree_ls = %s not recognized " % pDegree_ls
    elementQuadrature = CubeGaussQuadrature(nd, rotation_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(
        nd - 1, rotation_quad_order)
else:
    if pDegree_ls == 1:
        femSpaces = {0: C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ls == 2:
        femSpaces = {0: C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print "pDegree_ls = %s not recognized " % pDegree_ls
    elementQuadrature = SimplexGaussQuadrature(nd, rotation_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(
        nd - 1, rotation_quad_order)

subgridError = None

nLevels = ct.nLevels

subgridError = None  # HamiltonJacobi_ASGS_opt(coefficients, nd, lag=True)

massLumping = False


shockCapturing = NCLS.ShockCapturing(
    coefficients, nd, shockCapturingFactor=shockCapturingFactor_ls, lag=lag_shockCapturing_ls)

numericalFluxType = DoNothing  # StrongDirichlet


multilevelNonlinearSolver = Newton
# levelNonlinearSolver = Newton
#
#
# if ct.timeIntegration_ls == "be":
#     timeIntegration = BackwardEuler_cfl
#     stepController = Min_dt_controller
# elif ct.timeIntegration_ls == "vbdf":
#     timeIntegration = VBDF
#     stepController = Min_dt_cfl_controller
#     timeOrder = 2
# elif ct.timeIntegration_ls == "flcbdf":
#     timeIntegration = FLCBDF
#     stepController = FLCBDF_controller_sys
# elif ct.timeIntegration_ls == "rk":
#     if cDegree_ls == -1:
#         timeIntegration = LinearSSPRKPIintegration
#     else:
#         timeIntegration = LinearSSPRKintegration
#
#     timeIntegration = NCLS.RKEV
#     levelNonlinearSolverType = SSPRKNewton
#     stepController = Min_dt_RKcontroller
#     timeOrder = 1  # pDegree_ls + 1
#     nStagesTime = timeOrder
# else:
#     raise RuntimeError


# timeIntegration = NCLS.RKEV  # SSP33
# timeOrder = 1
# nStagesTime = timeOrder
# stepController = Min_dt_RKcontroller  # since substeps are substages

timeIntegration = BackwardEuler_cfl
# Serious error: it is not levelNonlinearSolverType
levelNonlinearSolver = ExplicitLumpedMassMatrix
stepController = Min_dt_cfl_controller


nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
nonlinearSmoother = None

fullNewtonFlag = False

tolFac = 0.0

nl_atol_res = atolLevelSet
l_atol_res = atolLevelSet

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py  # PETSc
    levelLinearSolver = KSP_petsc4py  # PETSc
    linear_solver_options_prefix = 'ncls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU

conservativeFlux = {}
maxNonlinearIts = 2
maxLineSearches = 1

#checkMass = True

if not applyCorrection and checkMass:
    auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryLS(
        "rotation2dnc" + `lRefinement`)]
