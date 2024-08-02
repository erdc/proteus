from proteus import *
from proteus.default_n import *
from la_periodicGauss_p import *

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
timeOrder=1
runCFL=0.33
triangleFlag = 1#if regular triangulatio then alternate diagonals

DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=41
nLevels = 1#4

if useVOS:
    subgridError = VOS3P.SubgridError(coefficients, nd)
    numericalFluxType = VOS3P.NumericalFlux
elif useNCLS:
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

if useVOS:
    shockCapturing = VOS3P.ShockCapturing(coefficients, nd, shockCapturingFactor=0.0, lag=True)
elif useNCLS:
    shockCapturing = NCLS.ShockCapturing(coefficients, nd, shockCapturingFactor=0.0, lag=True)
elif useVOF:
    shockCapturing = VOF.ShockCapturing(coefficients, nd, shockCapturingFactor=0.0, lag=True)
elif useHJ:
    shockCapturing = None#HamiltonJacobi_SC(coefficients,nd,lag=True)#None
else:
    shockCapturing = None#ResGrad_SC(coefficients,nd,lag=True)#None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

if useVOS:
    timeIntegration = VOS3P.RKEV
    timeOrder=3
#    levelNonlinearSolver = ExplicitLumpedMassMatrix
    levelNonlinearSolver = ExplicitConsistentMassMatrixForVOF
    fullNewtonFlag = True
    
nonlinearSmoother = NLGaussSeidel

tolFac = 0.0
nl_atol_res = 1.0e-8

linTolFac = 0.0
l_atol_res = 1.0e-8

maxNonlinearIts =1
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