from pyadh import *
from pyadh.default_n import *
from la_periodicGauss_2d_p import *


timeIntegration = FLCBDF
stepController = FLCBDF_controller
atol_u[0] = 1.0e-4
rtol_u[0] = 1.0e-4
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
runCFL=0.33
DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#nn=81
nLevels = 1#4

if useHJ:
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)
else:
    subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

numericalFluxType = DoNothing#None

if useHJ:
    shockCapturing = None#HamiltonJacobi_SC(coefficients,nd,lag=True)#None
else:
    shockCapturing = None#ResGrad_SC(coefficients,nd,lag=True)#None
    

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

maxNonlinearIts =3

matrix = SparseMatrix

multilevelLinearSolver = PETSc#LU

levelLinearSolver = PETSc#LU
#multilevelLinearSolver = LU

#levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
