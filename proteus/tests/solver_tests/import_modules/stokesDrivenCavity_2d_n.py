from proteus import *
from proteus.default_n import *
from proteus import defaults
defaults.reset_default_n()
try:
    from .stokesDrivenCavity_2d_p import *
except:
    from stokesDrivenCavity_2d_p import *

#######################################################

# context variables
nLevels = ct.nLevels
numeric_scheme = ct.numeric_scheme
useWeakBoundaryConditions = ct.useWeakBoundaryConditions
solveIteratively = ct.solveIteratively
solveInParallel = ct.solveInParallel

# TODO - this parallel flag needs work

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
    Q2 = LagrangeCubeFactory(2)
    femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis,
                 1:Q2,
                 2:Q2}
else:
    print('INVALID FINITE ELEMENT SELECTED')
#######################################################

######################################################
#           Mesh and Quadrature Options

if (numeric_scheme=="C0Q1C0Q1" or numeric_scheme=="THQuads"):
    # use a quadrilateral grid
    quad = True
    # numerical quadrature choice
    elementQuadrature = CubeGaussQuadrature(nd,2)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,1)
else:
    # use a simplex grid
    elementQuadrature = SimplexGaussQuadrature(nd,4)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
    # if unstructured would need triangleOptions flag to be set
    he = 0.05
    triangleOptions = "VApq30Dena%8.8f" % he
    
#####################################################

#no shock capturing
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
    numericalFluxType = Stokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #for true residual test
    linearSolverConvergenceTest = 'r-true'
    #to allow multiple models to set different ksp options
    #linear_solver_options_prefix = 'poisson_'
    linearSmoother = SimpleNavierStokes2D
else:
    if solveIteratively:
        multilevelLinearSolver = KSP_petsc4py
        levelLinearSolver = KSP_petsc4py
    else:
        multilevelLinearSolver = LU
        levelLinearSolver = LU
    if useWeakBoundaryConditions:
        numericalFluxType =  Stokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior 
    else:
        pass

#linear solver relative convergence test
linTolFac = 0.0
#linear solver absolute convergence test
l_atol_res = 1.0e-10

#linearSmoother=SimpleNavierStokes3D
#linearSmoother=NavierStokes3D_Qp
linearSmoother = petsc_LU
from petsc4py import PETSc
OptDB=PETSc.Options()
OptDB.clear()
OptDB.setValue("ksp_type","fgmres")
OptDB.setValue("ksp_atol",1e-20)
OptDB.setValue("ksp_atol",1e-12)
OptDB.setValue("pc_type","fieldsplit")
OptDB.setValue("pc_fieldsplit_type","schur")
OptDB.setValue("pc_fieldsplit_schur_fact_type","upper")
OptDB.setValue("fieldsplit_velocity_ksp_type","preonly")
OptDB.setValue("fieldsplit_velocity_pc_type","lu")
OptDB.setValue("fieldsplit_pressure_ksp_type","preonly")
