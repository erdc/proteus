from pyadh import *
from pyadh.default_n import *
from stokes_pome_ss_3d_p import *

timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

genMesh=True
he=0.25*radius
triangleOptions = "VpAq1.5ena%e" % (he**3 / 6.0)
#triangleOptions = "pAq1.5en"
nLevels = 1

subgridError = StokesASGS_velocity_pressure(coefficients,nd)
#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-10

matrix = SparseMatrix

multilevelLinearSolver = PETSc

levelLinearSolver = PETSc
#pick number of layers to use in overlap 
#nLayersOfOverlapForParallel = 0
#printLinearSolverInfo = True
#parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element

linearSmoother = GaussSeidel

linTolFac = 0.001

#conservativeFlux = None
#need exterior numerical flux for global conservation or parallel 
#numericalFluxType = Stokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior#None
numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior#None
#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior#None
