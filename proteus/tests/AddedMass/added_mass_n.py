from __future__ import absolute_import
from proteus.default_n import *
try:
    from . import added_mass_p as physics
except:
    import added_mass_p as physics
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.mprans import AddedMass
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions

elementQuadrature = ct.elementQuadrature
elementBoundaryQuadrature = ct.elementBoundaryQuadrature

triangleOptions = mesh.triangleOptions

femSpaces = {0:ct.basis}

stepController=StepControl.FixedStep

numericalFluxType =  AddedMass.NumericalFlux

matrix = LinearAlgebraTools.SparseMatrix


multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver      = LinearSolvers.KSP_petsc4py
parallelPartitioningType = mesh.parallelPartitioningType
nLayersOfOverlapForParallel = mesh.nLayersOfOverlapForParallel
nonlinearSmoother = NonlinearSolvers.AddedMassNewton
linearSmoother    = LinearSolvers.NavierStokesPressureCorrection

linear_solver_options_prefix = 'am_'

multilevelNonlinearSolver = NonlinearSolvers.AddedMassNewton
levelNonlinearSolver = NonlinearSolvers.AddedMassNewton

#linear solve rtolerance

linTolFac = 0.0
l_atol_res = ct.am_nl_atol_res
tolFac = 0.0
nl_atol_res = ct.am_nl_atol_res

nonlinearSolverConvergenceTest = 'rits'
levelNonlinearSolverConvergenceTest = 'rits'
linearSolverConvergenceTest             = 'r-true'
maxNonlinearIts =1
maxLineSearches=0
periodicDirichletConditions=None
conservativeFlux=None
auxiliaryVariables=[ct.system.ProtChAddedMass]
