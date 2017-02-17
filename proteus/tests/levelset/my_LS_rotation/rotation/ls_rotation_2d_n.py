from proteus import *
from proteus.default_n import *
from ls_rotation_2d_p import *
from rotation2D import *
nd = 2

multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton
fullNewtonFlag = fullNewton
updateJacobian = False

timeIntegration = NCLS.RKEV # SSP33 #mwf right now need timeIntegration to be SSP33 to run
#timeIntegration = BackwardEuler_cfl
stepController = Min_dt_cfl_controller

if timeIntegration_ncls == "SSP33": #mwf hack
    timeOrder = 3
    nStagesTime = 3
else:
    timeOrder = 1
    nStagesTime = 1

if useHex:
    if pDegree_ncls == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_ncls == 2:
        femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_ncls
    elementQuadrature = CubeGaussQuadrature(nd,rotation_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,rotation_quad_order)
else:
    if pDegree_ncls == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ncls == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_ncls
    elementQuadrature = SimplexGaussQuadrature(nd,rotation_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,rotation_quad_order)

nonlinearSmoother = None
subgridError = None
#subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=True)

numericalFluxType = NCLS.NumericalFlux
#numericalFluxType = DoNothing

shockCapturing = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_ncls,lag=lag_shockCapturing_ncls)

tolFac = 0.0

nl_atol_res = atolLevelSet
l_atol_res = atolLevelSet

maxNonlinearIts = 20
maxLineSearches = 0

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
if checkMass:
    auxiliaryVariables = [MassOverRegion()]

