from __future__ import absolute_import
from proteus import (StepControl,
                     TimeIntegration,
                     NonlinearSolvers,
                     LinearSolvers,
                     LinearAlgebraTools)
from proteus.default_n import *
from twp_navier_stokes_p import *

# MESH CONSTRUCTION #
triangleFlag = ct.triangleFlag if hasattr(ct,'triangleFlag') else None
nnx = ct.nnx if hasattr(ct,'nnx') else None
nny = ct.nny if hasattr(ct,'nny') else None
triangleOptions = ct.triangleOptions if hasattr(ct,'triangleOptions') else None

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
ns_shockCapturingFactor = ct.rans2p_parameters['ns_shockCapturingFactor']
ns_lag_shockCapturing = ct.rans2p_parameters['ns_lag_shockCapturing']
ns_lag_subgridError = ct.rans2p_parameters['ns_lag_subgridError']

#time stepping
timeDiscretization=ct.rans2p_parameters['timeDiscretization']
if timeDiscretization=='vbdf':
    timeIntegration = TimeIntegration.VBDF
    timeOrder = ct.rans2p_parameters['timeOrder']
    stepController  = StepControl.Min_dt_cfl_controller
else: #backward euler
    timeIntegration = TimeIntegration.BackwardEuler_cfl
    stepController  = StepControl.Min_dt_cfl_controller

elementQuadrature = ct.FESpace['elementQuadrature']
elementBoundaryQuadrature = ct.FESpace['elementBoundaryQuadrature']
if ct.nd==2:
    femSpaces = {0:ct.FESpace['basis'],
                 1:ct.FESpace['basis'],
                 2:ct.FESpace['basis']}
else:
    femSpaces = {0:ct.FESpace['basis'],
                 1:ct.FESpace['basis'],
                 2:ct.FESpace['basis'],
                 3:ct.FESpace['basis']}
    
massLumping       = False
numericalFluxType = RANS2P.NumericalFlux
subgridError = RANS2P.SubgridError(coefficients=coefficients,
                                   nd=nd,
                                   lag=ns_lag_subgridError,
                                   hFactor=ct.FESpace['hFactor'])
shockCapturing = RANS2P.ShockCapturing(coefficients=coefficients,
                                       nd=nd,
                                       shockCapturingFactor=ns_shockCapturingFactor,
                                       lag=ns_lag_shockCapturing)

fullNewtonFlag = True
multilevelNonlinearSolver = NonlinearSolvers.Newton
levelNonlinearSolver      = NonlinearSolvers.Newton

nonlinearSmoother = None
linearSmoother    = LinearSolvers.SimpleNavierStokes2D
matrix = LinearAlgebraTools.SparseMatrix

multilevelLinearSolver = LinearSolvers.KSP_petsc4py
levelLinearSolver      = LinearSolvers.KSP_petsc4py
if ct.opts.useSuperlu:
    multilevelLinearSolver = LinearSolvers.LU
    levelLinearSolver      = LinearSolvers.LU

linear_solver_options_prefix = 'rans2p_'
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

# TOLERANCS #
ns_nl_atol_res = max(1.0e-8, 0.001 * ct.he ** 2)
nl_atol_res = ns_nl_atol_res
tolFac = 0.0
linTolFac = 0.01
l_atol_res = 0.01*ns_nl_atol_res

useEisenstatWalker = False
maxNonlinearIts = 50
maxLineSearches = 0
conservativeFlux = {0:'pwl-bdm-opt'}
