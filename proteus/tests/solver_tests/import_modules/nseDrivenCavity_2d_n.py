from proteus import *
from proteus.default_n import *
from proteus import defaults
defaults.reset_default_n()
try:
    from .nseDrivenCavity_2d_p import *
except:
    from nseDrivenCavity_2d_p import *
import petsc4py

#######################################################

# context variables
nLevels = ct.nLevels
numeric_scheme = ct.numeric_scheme
useWeakBoundaryConditions = ct.useWeakBoundaryConditions
solveIteratively = ct.solveIteratively
solveInParallel = ct.solveInParallel
schur_solver = ct.schur_solver

#######################################################

#steady-state so no time integration
stepController = FixedStep
timeIntegration = NoIntegration
#number of output timesteps

import numpy as  np
tnList =  [0.0,1.0,1.7,2.2,2.7,3.2,3.5,3.7]#,1.5,2.0,2.5,3.0]

shockCapturing = None

######################################################
#            specify the numerical scheme
#no shock capturing

if numeric_scheme=="TH":
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 2:C0_AffineQuadraticOnSimplexWithNodalBasis}
elif numeric_scheme=="C0P1C0P1":
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis}
    subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
elif numeric_scheme=="C0Q1C0Q1":
    femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis,
                1:C0_AffineLinearOnCubeWithNodalBasis,
                2:C0_AffineLinearOnCubeWithNodalBasis}
    subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
elif numeric_scheme=="THQuads":
    femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis,
                 1:Q2,
                 2:Q2}
else:
    print('INVALID FINITE ELEMENT SELECTED')
#######################################################
#              Mesh and Quadrature Options

if (numeric_scheme=="C0Q1C0Q1" or numeric_scheme=="THQuads"):
    #use a quadrilateral grid
    quad = True
    # numerical quadrature choice
    elementQuadrature = CubeGaussQuadrature(nd,3)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
else:
    #use a simplex grid
    elementQuadrature = SimplexGaussQuadrature(nd,4)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
    #if unstructured would need triangleOptions flag to be set
    triangleOptions = "VApq30Dena%8.8f" % he


#nonlinear solver choices
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
#linear problem so force 1 iteration allowed
maxNonlinearIts = 25
maxLineSearches = 0
fullNewtonFlag = True
#absolute nonlinear solver residual tolerance
#nl_atol_res = 1.0e-8
nl_atol_res = 1.0e-5
nl_rtol_res = 1.0e-20
#nl_atol_res = 1.0-5
#nl_rtol_res = 1.0e-6
#relative nonlinear solver convergence tolerance as a function of h
#(i.e., tighten relative convergence test as we refine)
tolFac = 0.

#matrix type
matrix = SparseMatrix

if solveInParallel:
    multilevelLinearSolver = KSP_petsc4py
    #for petsc do things like
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = KSP_petsc4py#
    #for petsc do things like
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    #levelLinearSolver = PETSc#
    #pick number of layers to use in overlap
    nLayersOfOverlapForParallel = 0
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    #have to have a numerical flux in parallel
    numericalFluxType = NavierStokes_AdvectionDiagonalUpwinde_Diffusion_IIPG_exterior
    #for true residual test
    linearSolverConvergenceTest = 'r-true'
    #to allow multiple models to set different ksp options
    #linear_solver_options_prefix = 'poisson_'
    linearSmoother = None
else:
    if solveIteratively:
        linearSolverConvergenceTest = 'r-true'
        multilevelLinearSolver = KSP_petsc4py
        levelLinearSolver = KSP_petsc4py
    else:
        multilevelLinearSolver = LU
        levelLinearSolver = LU
    if useWeakBoundaryConditions:
        numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior
        #numericalFluxType = RANS2P.NumericalFlux
        numericalFluxType.useStrongDirichletConstraints = True
    else:
        pass
    
linTolFac = 0.001
#linear solver absolute convergence test
l_atol_res = 1.0e-8

if schur_solver == 'Qp':
    linearSmoother = Schur_Qp
elif schur_solver == 'LSC':
    linearSmoother=Schur_LSC
elif schur_solver == 'selfp':
    linearSmoother=SimpleNavierStokes3D
elif schur_solver == 'petsc_LU':
    linearSmoother=petsc_LU
else:
    raise Exception('invalid solver type')
