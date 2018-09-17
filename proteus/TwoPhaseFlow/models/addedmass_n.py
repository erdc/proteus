from __future__ import absolute_import
from proteus.default_n import *
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.mprans import AddedMass
import addedmass_p as physics

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
triangleFlag = ct.triangleFlag if hasattr(ct,'triangleFlag') else None
nnx = ct.nnx if hasattr(ct,'nnx') else None
nny = ct.nny if hasattr(ct,'nny') else None
nnz = ct.nnz if hasattr(ct,'nnz') else None
triangleOptions = ct.triangleOptions if hasattr(ct,'triangleOptions') and ct.triangleOptions != 'q30DenA' else mesh.triangleOptions
if hasattr(ct, 'nLevels'):
    nLevels = ct.nLevels
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
restrictFineSolutionToAllMeshes = mesh.restrictFineSolutionToAllMeshes

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
stepController = StepControl.FixedStep

# ******************************************* #
# ********** FINITE ELEMENT SPACES ********** #
# ******************************************* #
elementQuadrature = ct.FESpace['elementQuadrature']
elementBoundaryQuadrature = ct.FESpace['elementBoundaryQuadrature']
femSpaces = {0: ct.FESpace['lsBasis']}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver = NonlinearSolvers.AddedMassNewton
levelNonlinearSolver = NonlinearSolvers.AddedMassNewton
nonlinearSmoother = NonlinearSolvers.AddedMassNewton
#
nonlinearSolverConvergenceTest = 'r'
levelNonlinearSolverConvergenceTest = 'r'

# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
numericalFluxType = AddedMass.NumericalFlux
conservativeFlux = None

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = LinearAlgebraTools.SparseMatrix
linearSmoother = LinearSolvers.NavierStokesPressureCorrection
multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver = LinearSolvers.KSP_petsc4py
if ct.opts.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver = LinearSolvers.LU
#
linear_solver_options_prefix = 'am_'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = max(ct.params.minTol, ct.params.addedmass['tolFac']*ct.he**2)
linTolFac = 0.
l_atol_res = nl_atol_res
#
tolFac = 0.
maxNonlinearIts = 1
maxLineSearches = 0


auxiliaryVariables=[ct.system.ProtChAddedMass]
