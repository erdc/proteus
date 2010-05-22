from pyadh import *
from pyadh.default_n import *
from stokes_cavity_ss_2d_p import *

timeIntegration = NoIntegration
timeIntegrator = SteadyStateIntegrator
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 51
nLevels = 1

subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#PETSc eg #-pc_type lu -pc_factor_mat_solver_package
multilevelLinearSolver = PETSc#LU
levelLinearSolver = PETSc#LU
#multilevelLinearSolver = LU
#levelLinearSolver = LU

# multilevelLinearSolver = StarILU
# multilevelLinearSolver = NI
# levelLinearSolver = StarILU
# linearSmoother = StarILU

# multilevelLinearSolver = LU
# levelLinearSolver = LU
# linearSmoother = StarILU

printLinearSolverInfo = True
linTolFac = 0.0

numericalFluxType = Stokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation
conservativeFlux = None
