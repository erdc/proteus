from proteus import *
from proteus.default_n import *
from la_periodicGauss_p import *

timeIntegration = VBDF
stepController = Min_dt_controller
timeOrder=2
runCFL=0.33

DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=11
nLevels = 1#4

if useNCLS:
    subgridError = NCLS.SubgridError(coefficients, nd)
    numericalFluxType = NCLS.NumericalFlux
elif useVOF:
    subgridError = VOF.SubgridError(coefficients, nd)
    numericalFluxType = VOF.NumericalFlux
elif useHJ:
    subgridError = HamiltonJacobi_ASGS(coefficients, nd, lag=False)
    numericalFluxType = DoNothing#None
else:
    subgridError = AdvectionDiffusionReaction_ASGS(coefficients, nd, lag=False)
    numericalFluxType = DoNothing#None

if useNCLS:
    shockCapturing = NCLS.ShockCapturing(coefficients, nd, shockCapturingFactor=0.0, lag=True)
elif useVOF:
    shockCapturing = VOF.ShockCapturing(coefficients, nd, shockCapturingFactor=0.0, lag=True)
elif useHJ:
    shockCapturing = None#HamiltonJacobi_SC(coefficients,nd,lag=True)#None
else:
    shockCapturing = None#ResGrad_SC(coefficients,nd,lag=True)#None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0
nl_atol_res = 1.0e-8

linTolFac = 0.0
l_atol_res = 1.0e-8

maxNonlinearIts =3
maxLineSearches = 0

matrix = SparseMatrix

multilevelLinearSolver = LU
#
levelLinearSolver = LU

#multilevelLinearSolver = KSP_petsc4py

#levelLinearSolver = KSP_petsc4py

linearSmoother = None

conservativeFlux = None

checkMass = True

parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0
