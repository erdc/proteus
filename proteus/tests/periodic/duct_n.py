from proteus import *
from proteus.default_n import *
from duct_p import *

#time stepping
runCFL = 0.33
timeIntegration = VBDF
timeOrder = 2
stepController  = StepControl.Min_dt_cfl_controller
systemStepExact = False

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#nnx=81
#nny=21
#nnz=21
#nnx=41
#nny=11
#nnz=11
nnx=21
nny=6
nnz=6
#nLevels = 1

numericalFluxType = RANS2P.NumericalFlux

parallelPeriodic=True

subgridError = RANS2P.SubgridError(coefficients=coefficients,
                                  nd=nd,
                                  lag=True,
                                  hFactor=1.0)

shockCapturing = RANS2P.ShockCapturing(coefficients=coefficients,
                                      nd=nd,
                                      shockCapturingFactor=0.0,
                                      lag=True)

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True
maxNonlinearIts = 50
maxLineSearches = 0

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = KSP_petsc4py
levelLinearSolver = KSP_petsc4py
linearSmoother = SimpleNavierStokes3D
#multilevelLinearSolver = LU
#levelLinearSolver = LU
linearSmoother = None

linear_solver_options_prefix = 'rans2p_'

linTolFac = 0.001

conservativeFlux = None

parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
