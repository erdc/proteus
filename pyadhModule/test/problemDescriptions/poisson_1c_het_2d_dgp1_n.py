from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_2d_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = dict((i,DG_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 11
nLevels = 2

subgridError = None

shockCapturing = None

numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG
#numericalFluxType = Diffusion_LDG

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-6

nl_atol_res = 1.0e-12

matrix = SparseMatrix

multilevelLinearSolver = NI
multilevelLinearSolver = LU

levelLinearSolver = LU
levelLinearSolver = LU

linearSmoother = StarILU#GaussSeidel#Jacobi#StarILU

linTolFac = 1.0e-8

cfluxtag  = 'dg-bdm' #'dg-point-eval','dg','dg-bdm'
conservativeFlux = dict((i,cfluxtag) for i in range(nc))
multigridCycles = 3

preSmooths = 3

postSmooths = 3

archiveFlag = ArchiveFlags.EVERY_USER_STEP
