from proteus import *
from redist_p import *
from tank import *

nl_atol_res = rd_nl_atol_res
tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.1*rd_nl_atol_res

if redist_Newton:
    timeIntegration = NoIntegration
    stepController = Newton_controller
    maxNonlinearIts = 25
    maxLineSearches = 0
    nonlinearSolverConvergenceTest = 'r'
    levelNonlinearSolverConvergenceTest = 'r'
    linearSolverConvergenceTest = 'r-true'
    useEisenstatWalker = True
else:
    timeIntegration = BackwardEuler_cfl
    stepController = RDLS.PsiTC
    runCFL=1.0
    psitc['nStepsForce']=3
    psitc['nStepsMax']=25
    psitc['reduceRatio']=2.0
    psitc['startRatio']=1.0
    rtol_res[0] = 0.0
    atol_res[0] = rd_nl_atol_res
    useEisenstatWalker = False
    maxNonlinearIts = 1
    maxLineSearches = 0
    nonlinearSolverConvergenceTest = 'rits'
    levelNonlinearSolverConvergenceTest = 'rits'
    linearSolverConvergenceTest = 'r-true'

femSpaces = {0:basis}

massLumping       = False
numericalFluxType = DoNothing
conservativeFlux  = None
subgridError      = RDLS.SubgridError(coefficients,nd)
shockCapturing    = RDLS.ShockCapturing(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=rd_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver  = Newton
levelNonlinearSolver       = Newton

nonlinearSmoother = NLGaussSeidel
linearSmoother    = None

matrix = SparseMatrix

if useOldPETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver      = PETSc
else:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver      = KSP_petsc4py

if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU

linear_solver_options_prefix = 'rdls_'

#auxiliaryVariables=[lineGauges_phi]
