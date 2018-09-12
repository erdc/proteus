from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from rans3p_p import *

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
triangleFlag = ct.triangleFlag if hasattr(ct,'triangleFlag') else None
nnx = ct.nnx if hasattr(ct,'nnx') else None
nny = ct.nny if hasattr(ct,'nny') else None
triangleOptions = ct.triangleOptions if hasattr(ct,'triangleOptions') and ct.triangleOptions != 'q30DenA' else mesh.triangleOptions

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
ns_shockCapturingFactor = ct.rans3p_parameters['ns_shockCapturingFactor']
ns_lag_shockCapturing = ct.rans3p_parameters['ns_lag_shockCapturing']
ns_lag_subgridError = ct.rans3p_parameters['ns_lag_subgridError']

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeDiscretization=ct.rans3p_parameters['timeDiscretization']
if timeDiscretization=='vbdf':
    timeIntegration = VBDF
    timeOrder=2
    stepController  = Min_dt_cfl_controller
else: #backward euler
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_cfl_controller
runCFL=ct.opts.cfl

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = ct.FESpace['elementQuadrature']
elementBoundaryQuadrature = ct.FESpace['elementBoundaryQuadrature']
if nd==2:
    femSpaces = {0:ct.FESpace['velBasis'],
	         1:ct.FESpace['velBasis']}
else:
    femSpaces = {0:ct.FESpace['velBasis'],
	         1:ct.FESpace['velBasis'],
                 2:ct.FESpace['velBasis']}

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
                                    hFactor=ct.FESpace['hFactor'])
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
if ct.opts.useSuperlu:
    multilevelLinearSolver = LU
    levelLinearSolver      = LU
#
linear_solver_options_prefix = 'rans3p_'
linearSolverConvergenceTest             = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
ns_nl_atol_res = max(1.0e-10, 0.01 * ct.he ** 2)
nl_atol_res = ns_nl_atol_res
tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ns_nl_atol_res

useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
