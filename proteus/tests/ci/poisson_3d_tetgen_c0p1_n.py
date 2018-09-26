from __future__ import absolute_import
from __future__ import division
from past.utils import old_div
from proteus import *
from proteus.default_n import *
try:
    from .poisson_3d_tetgen_p import *
except:
    from poisson_3d_tetgen_p import *

#steady-state so no time integration
timeIntegration = NoIntegration
#number of output timesteps
nDTout = 1

#finite element spaces
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
#numerical quadrature choices
elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

triangleOptions="VApq1.35q12feena%e" % (old_div((he**3),6.0),)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))

#number of levels in mesh
nLevels = 1

#no stabilization or shock capturing
subgridError = None

shockCapturing = None

#nonlinear solver choices
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
#linear problem so force 1 iteration allowed
maxNonlinearIts = 2
maxLineSearches = 1
fullNewtonFlag = True
#absolute nonlinear solver residual tolerance
nl_atol_res = 1.0e-8
#relative nonlinear solver convergence tolerance as a function of h
#(i.e., tighten relative convergence test as we refine)
tolFac = 0.0

#matrix type
matrix = SparseMatrix

#convenience flag
parallel = True

if parallel:
    multilevelLinearSolver = KSP_petsc4py
    #for petsc do things lie
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
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
    #for true residual test
    linearSolverConvergenceTest = 'r-true'
    #to allow multiple models to set different ksp options
    #linear_solver_options_prefix = 'poisson_'
    linearSmoother = None
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

#linear solver relative convergence test
linTolFac = 0.0
#linear solver absolute convergence test
l_atol_res = 1.0e-10

conservativeFlux =  None
cfluxtag = None
