from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
from pressureincrement_p import *

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
triangleFlag = ct.triangleFlag if hasattr(ct,'triangleFlag') else None
nnx = ct.nnx if hasattr(ct,'nnx') else None
nny = ct.nny if hasattr(ct,'nny') else None
triangleOptions = ct.triangleOptions if hasattr(ct,'triangleOptions') and ct.triangleOptions != 'q30DenA' else mesh.triangleOptions

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
stepController=FixedStep

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = ct.FESpace['elementQuadrature']
elementBoundaryQuadrature = ct.FESpace['elementBoundaryQuadrature']
femSpaces = {0:ct.FESpace['pBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver = NonlinearSolvers.Newton
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
numericalFluxType = PresInc.NumericalFlux
conservativeFlux=None

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver = LinearSolvers.KSP_petsc4py

linearSmoother = None
if ct.opts.useSuperlu:
    linearSmoother    = None
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'phi_'
linearSolverConvergenceTest             = 'r-true'
maxLineSearches=0

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
pressure_nl_atol_res = max(1.0e-10, 0.01 * ct.he ** 2)
nl_atol_res = pressure_nl_atol_res
tolFac = 0.0
linTolFac = 0.0
l_atol_res = 0.01*pressure_nl_atol_res
maxNonlinearIts = 50
