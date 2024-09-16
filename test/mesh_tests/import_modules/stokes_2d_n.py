from proteus import *
from proteus.default_n import *
try:
    from .stokes_2d_p import *
except:
    from stokes_2d_p import *

#######################################################

# context variables
nLevels = 2
numeric_scheme = "THQuads"
useWeakBoundaryConditions = False
solveIteratively = False
solveInParallel = False

#######################################################

#steady-state so no time integration
timeIntegration = NoIntegration
#number of output timesteps
nDTout = 1

######################################################
#            specify the numerical scheme
if numeric_scheme=="TH":
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineQuadraticOnSimplexWithNodalBasis,
                 2:C0_AffineQuadraticOnSimplexWithNodalBasis}
elif numeric_scheme=="C0P1C0P1":
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
                 1:C0_AffineLinearOnSimplexWithNodalBasis,
                 2:C0_AffineLinearOnSimplexWithNodalBasis}
    subgridError = StokesASGS_velocity_pressure(coefficients,nd)
elif numeric_scheme=="C0Q1C0Q1":
    femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis,
                 1:C0_AffineLinearOnCubeWithNodalBasis,
                 2:C0_AffineLinearOnCubeWithNodalBasis}
    subgridError = StokesASGS_velocity_pressure(coefficients,nd)
elif numeric_scheme=="THQuads":
    femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis,
                 1:Q2,
                 2:Q2}
else:
    print('INVALID FINITE ELEMENT SELECTED')
#######################################################
#              Mesh and Quadrature Options

he = 0.005
if (numeric_scheme=="C0Q1C0Q1" or numeric_scheme=="THQuads"):
    #use a quadrilateral grid
    quad = True
    # numerical quadrature choice
    elementQuadrature = CubeGaussQuadrature(nd,4)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    #nnx = nny  = int(ceil(L[0]/he) + 1)
else:
    #use a simplex grid
    elementQuadrature = SimplexGaussQuadrature(nd,4)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
    #if unstructured would need triangleOptions flag to be set
    triangleOptions = "VApq30Dena%8.8f" % he
#no stabilization or shock capturing

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

if solveInParallel:
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
    numericalFluxType = Stokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #for true residual test
    linearSolverConvergenceTest = 'r-true'
    #to allow multiple models to set different ksp options
    #linear_solver_options_prefix = 'poisson_'
    linearSmoother = None
else:
#    numericalFluxType = Diffusion_IIPG_exterior
    if solveIteratively:
        multilevelLinearSolver = KSP_petsc4py
        levelLinearSolver = KSP_petsc4py
    else:
        multilevelLinearSolver = LU
        levelLinearSolver = LU
    if useWeakBoundaryConditions:
        numericalFluxType = Stokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior
#        numericalFluxType = Diffusion_IIPG_exterior
    else:
        pass

#linear solver relative convergence test
linTolFac = 0.0
linearSmoother = Schur_Qp
#linear solver absolute convergence test
l_atol_res = 1.0e-10

#conservativeFlux =  {0:'pwl'}
conservativeFlux =  None
cfluxtag = None

hex=False
