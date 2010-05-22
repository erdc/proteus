from pyadh import *
from pyadh.default_n import *
#make sure to change name of imported p file when you copy this file
from navier_stokes_poiseuille_ss_3d_p import *

timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

nnx=3
nny=3
nnz=3
#number of refinement  levels, might want to decrease for 3D

nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts =10#1#100

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-7

matrix = SparseMatrix

multilevelLinearSolver = LU
levelLinearSolver = LU
multilevelLinearSolver = PETSc
levelLinearSolver = PETSc

linTolFac = 0.001

# conservativeFlux = {0:'pwl',
#                     1:'pwl',
#                     2:'pwl'}

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
