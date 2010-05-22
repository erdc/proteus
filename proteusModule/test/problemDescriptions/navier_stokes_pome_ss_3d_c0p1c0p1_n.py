from pyadh import *
from pyadh.default_n import *
from navier_stokes_pome_ss_3d_p import *

timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

triangleOptions = "Aen"
nn=11#3
nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True) #None

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU
levelLinearSolver = LU
multilevelLinearSolver = PETSc
levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior

parallelPartitioningType = pyadh.MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0 
