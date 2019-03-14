from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from rans3p_p import *

# *********************************************** #
# ********** Read from myTpFlowProblem ********** #
# *********************************************** #
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

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
ns_shockCapturingFactor = rans3p_parameters['ns_shockCapturingFactor']
ns_lag_shockCapturing = rans3p_parameters['ns_lag_shockCapturing']
ns_lag_subgridError = rans3p_parameters['ns_lag_subgridError']

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeDiscretization=rans3p_parameters['timeDiscretization']
if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
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
subgridError = RANS3PF.SubgridError(coefficients=coefficients,
                                    nd=nd,
                                    lag=ns_lag_subgridError,
                                    hFactor=FESpace['hFactor'])
shockCapturing = RANS3PF.ShockCapturing(coefficients=coefficients,
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
