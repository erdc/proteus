from proteus import *
from proteus.default_n import *
from md_flow_p import *

parallel = False
timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
maxNonlinearIts = 1 #still linear in pressure/head

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py
    #for petsc do things lie
    #"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10" or
    #-pc_type lu -pc_factor_mat_solver_package
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    levelLinearSolver = KSP_petsc4py
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG_exterior
    #to allow multiple models to set different ksp options
    linear_solver_options_prefix = 'flow_'
    linearSmoother = None
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    linearSolverConvergenceTest= 'r'

linearSmoother = StarILU

linTolFac = 0.0
l_atol_res = 1.0e-10

cfluxtag  = 'p1-nc'
conservativeFlux =  {0:cfluxtag}
