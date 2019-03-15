from __future__ import absolute_import
<<<<<<< HEAD
from proteus import *
from proteus.default_n import *
from rans2p_p import *
=======
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.mprans import RANS2P
import rans2p_p as physics
>>>>>>> TwoPhaseFlow

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
<<<<<<< HEAD
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
=======
ct = physics.ct
myTpFlowProblem = physics.myTpFlowProblem
nd = myTpFlowProblem.nd
cfl = myTpFlowProblem.cfl
FESpace = myTpFlowProblem.FESpace
useSuperlu = myTpFlowProblem.useSuperlu
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
myparams = params.Models.rans2p
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
>>>>>>> TwoPhaseFlow

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
<<<<<<< HEAD
ns_shockCapturingFactor = rans2p_parameters['ns_shockCapturingFactor']
ns_lag_shockCapturing = rans2p_parameters['ns_lag_shockCapturing']
ns_lag_subgridError = rans2p_parameters['ns_lag_subgridError']
=======
ns_shockCapturingFactor = myparams.ns_shockCapturingFactor
ns_lag_shockCapturing = myparams.ns_lag_shockCapturing
ns_lag_subgridError = myparams.ns_lag_subgridError
>>>>>>> TwoPhaseFlow

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
<<<<<<< HEAD
timeDiscretization=rans2p_parameters['timeDiscretization']
if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder = rans2p_parameters['timeOrder']
=======
timeDiscretization = myparams.timeDiscretization
if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder = myparams.timeOrder
>>>>>>> TwoPhaseFlow
    stepController  = Min_dt_cfl_controller
else: #backward euler
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller
runCFL=cfl

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = FESpace['elementQuadrature']
elementBoundaryQuadrature = FESpace['elementBoundaryQuadrature']
if nd==2:
    femSpaces = {0: FESpace['pBasis'],
                 1: FESpace['velBasis'],
                 2: FESpace['velBasis']}
else:
    femSpaces = {0: FESpace['pBasis'],
                 1: FESpace['velBasis'],
                 2: FESpace['velBasis'],
                 3: FESpace['velBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
fullNewtonFlag = True
<<<<<<< HEAD
multilevelNonlinearSolver = Newton
levelNonlinearSolver      = Newton
=======
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton
>>>>>>> TwoPhaseFlow
nonlinearSmoother = None
nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'r'

# ******************************************************** #
# ********** NUMERICAL FLUXES AND STABILIZATION ********** #
# ******************************************************** #
numericalFluxType = RANS2P.NumericalFlux
<<<<<<< HEAD
conservativeFlux = None #{0:'pwl-bdm-opt'}
subgridError = RANS2P.SubgridError(coefficients=coefficients,
                                   nd=nd,
                                   lag=ns_lag_subgridError,
                                   hFactor=FESpace['hFactor'])
shockCapturing = RANS2P.ShockCapturing(coefficients=coefficients,
=======
# conservativeFlux = {0:'pwl-bdm-opt'}
conservativeFlux = None
subgridError = RANS2P.SubgridError(coefficients=physics.coefficients,
                                   nd=nd,
                                   lag=ns_lag_subgridError,
                                   hFactor=FESpace['hFactor'])
shockCapturing = RANS2P.ShockCapturing(coefficients=physics.coefficients,
>>>>>>> TwoPhaseFlow
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
linear_solver_options_prefix = 'rans2p_'
linearSolverConvergenceTest             = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
<<<<<<< HEAD
ns_nl_atol_res = max(1.0e-8, 0.01 * he ** 2)
=======
ns_nl_atol_res = max(myparams.minTol, myparams.tolFac*he**2)
>>>>>>> TwoPhaseFlow
nl_atol_res = ns_nl_atol_res
tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ns_nl_atol_res

useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
<<<<<<< HEAD
=======

auxiliaryVariables = myparams.auxiliaryVariables
>>>>>>> TwoPhaseFlow
