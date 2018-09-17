from __future__ import absolute_import
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.mprans import VOF
import vof_p as physics

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
ct = physics.ct
myTpFlowProblem = physics.myTpFlowProblem
nd = myTpFlowProblem.nd
params = myTpFlowProblem.Parameters
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
he = myTpFlowProblem.he
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

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
runCFL = ct.opts.cfl

# ******************************************* #
# ********** FINITE ELEMENT SPACES ********** #
# ******************************************* #
elementQuadrature = ct.FESpace['elementQuadrature']
elementBoundaryQuadrature = ct.FESpace['elementBoundaryQuadrature']
femSpaces = {0: ct.FESpace['lsBasis']}

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
numericalFluxType = VOF.NumericalFlux
conservativeFlux = None
subgridError = VOF.SubgridError(coefficients=physics.coefficients,
                                nd=nd)
shockCapturing = VOF.ShockCapturing(coefficients=physics.coefficients,
                                    nd=nd,
                                    shockCapturingFactor=params.vof['shockCapturingFactor'],
                                    lag=params.vof['lag_shockCapturing'])

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
linearSmoother = None
multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver = LinearSolvers.KSP_petsc4py
if ct.opts.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'vof_'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = max(params.minTol, params.vof['tolFac']*ct.he**2)
linTolFac = 0.001
l_atol_res = 0.001*nl_atol_res
#
useEisenstatWalker = False#True
tolFac = 0.
maxNonlinearIts = 50
maxLineSearches = 0


auxiliaryVariables = ct.domain.auxiliaryVariables['vof']
