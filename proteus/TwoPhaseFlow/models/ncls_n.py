from __future__ import absolute_import
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.mprans import NCLS
import ncls_p as physics

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
ct = physics.ct
myTpFlowProblem = physics.myTpFlowProblem
nd = myTpFlowProblem.nd
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
he = myTpFlowProblem.he
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
triangleFlag = myTpFlowProblem.triangleFlag
nnx = myTpFlowProblem.nnx
nny = myTpFlowProblem.nny
nnz = myTpFlowProblem.nnz
triangleOptions = domain.MeshOptions.triangleOptions

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeIntegration = TimeIntegration.BackwardEuler_cfl
stepController  = StepControl.Min_dt_cfl_controller
runCFL = cfl

# ******************************************* #
# ********** FINITE ELEMENT SPACES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
femSpaces = {0: FESpace['lsBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton
fullNewtonFlag = True
nonlinearSmoother = None
#
levelNonlinearSolverConvergenceTest = 'r'

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
massLumping = False
numericalFluxType = NCLS.NumericalFlux
conservativeFlux = None
subgridError = NCLS.SubgridError(coefficients=physics.coefficients,
                                nd=nd)
shockCapturing = NCLS.ShockCapturing(coefficients=physics.coefficients,
                                     nd=nd,
                                     shockCapturingFactor=mparams.ncls['shockCapturingFactor'],
                                     lag=mparams.ncls['lag_shockCapturing'])

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
linearSmoother = None
multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver = LinearSolvers.KSP_petsc4py
if useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'ncls_'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = max(mparams.ncls.minTol, mparams.ncls['tolFac']*he**2)
linTolFac = 0.001
l_atol_res = 0.001*nl_atol_res
#
useEisenstatWalker = False#True
tolFac = 0.
maxNonlinearIts = 50
maxLineSearches = 0

auxiliaryVariables = domain.auxiliaryVariables['ls']
