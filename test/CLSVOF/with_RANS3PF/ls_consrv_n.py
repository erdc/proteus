from proteus import *
from .multiphase import *
from .ls_consrv_p import *

timeIntegrator  = ForwardIntegrator
timeIntegration = NoIntegration

femSpaces = {0:pbasis}

subgridError      = None
massLumping       = False
numericalFluxType = DoNothing
conservativeFlux  = None
shockCapturing    = None

fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton

nonlinearSmoother = None
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

linear_solver_options_prefix = 'mcorr_'
nonlinearSolverConvergenceTest  = 'rits'
levelNonlinearSolverConvergenceTest  = 'rits'
linearSolverConvergenceTest  = 'r-true'

tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*mcorr_nl_atol_res
nl_atol_res = mcorr_nl_atol_res
useEisenstatWalker = False

maxNonlinearIts = 50
maxLineSearches = 0
