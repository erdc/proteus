from proteus import *
from proteus.default_n import *
from .sp_gw_p import *
reload(proteus.default_n)

#fixed time step
timeIntegration = BackwardEuler
DT = float(T/nDTout)
stepController = FixedStep

#Piecewise-Linears on triangles
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
maxNonlinearIts = 1 #linear in pressure/head

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
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    linearSolverConvergenceTest= 'r'

#
linTolFac = 0.0
l_atol_res = 1.0e-10

#for post-processing velocities to get locally conservative approximations
cfluxtag  = 'pwl'
conservativeFlux =  {0:cfluxtag}

