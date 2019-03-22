from __future__ import absolute_import
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.mprans import RANS3PF
import rans3p_p as physics

forceTerms = physics.forceTerms

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
ct = physics.ct
myTpFlowProblem = physics.myTpFlowProblem
nd = myTpFlowProblem.nd
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
myparams = mparams.rans3p
pparams = params.physical # physical parameters
meshparams = params.mesh

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
he = meshparams.he
triangleFlag = meshparams.triangleFlag
nnx = meshparams.nnx
nny = meshparams.nny
nnz = meshparams.nnz
triangleOptions = meshparams.triangleOptions
parallelPartitioningType = meshparams.parallelPartitioningType
nLayersOfOverlapForParallel = meshparams.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = meshparams.restrictFineSolutionToAllMeshes

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
ns_shockCapturingFactor = myparams.ns_shockCapturingFactor
ns_lag_shockCapturing = myparams.ns_lag_shockCapturing
ns_lag_subgridError = myparams.ns_lag_subgridError

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeDiscretization = myparams.timeDiscretization
if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
    stepController  = Min_dt_cfl_controller
else: #backward euler
    timeIntegration = TimeIntegration.BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller
runCFL=cfl

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
if nd==2:
    femSpaces = {0:FESpace['velBasis'],
	         1:FESpace['velBasis']}
else:
    femSpaces = {0:FESpace['velBasis'],
	         1:FESpace['velBasis'],
                 2:FESpace['velBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
fullNewtonFlag = True
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton
nonlinearSmoother = None
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'

# ******************************************************** #
# ********** NUMERICAL FLUXES AND STABILIZATION ********** #
# ******************************************************** #
numericalFluxType = RANS3PF.NumericalFlux
conservativeFlux  = None
subgridError = RANS3PF.SubgridError(coefficients=physics.coefficients,
                                    nd=nd,
                                    lag=ns_lag_subgridError,
                                    hFactor=FESpace['hFactor'])
shockCapturing = RANS3PF.ShockCapturing(coefficients=physics.coefficients,
                                        nd=nd,
                                        shockCapturingFactor=ns_shockCapturingFactor,
                                        lag=ns_lag_shockCapturing)

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
massLumping = False
matrix = SparseMatrix
linearSmoother = None
# Linear solver 
multilevelLinearSolver = KSP_petsc4py
levelLinearSolver      = KSP_petsc4py
if useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU
#
linear_solver_options_prefix = 'rans3p_'
linearSolverConvergenceTest             = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
ns_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
nl_atol_res = ns_nl_atol_res
tolFac = 0.0
linTolFac = 0.0
l_atol_res = 0.1*ns_nl_atol_res

useEisenstatWalker = False
maxNonlinearIts = 1 # This is a linear problem
maxLineSearches = 0

auxiliaryVariables = myparams.auxiliaryVariables
