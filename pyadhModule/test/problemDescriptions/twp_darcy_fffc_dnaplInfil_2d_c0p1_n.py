from pyadh import *
from pyadh.default_n import *
from twp_darcy_fffc_dnaplInfil_2d_p import *

#timeIntegration = BackwardEuler
#timeIntegrator = ForwardIntegrator
#runCFL=0.1
#if useForsyth:
#    DT = 1.0e2#1.0e-2
#else:
#    DT = 10.0
runCFL=None
timeIntegration = FLCBDF
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-4
atol_u[0] = 1.0e-4
rtol_u[1] = 1.0e-4
atol_u[1] = 1.0e-4
DT = None
nDTout = 50#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=41
#nLevels = 1

#nn=3
#nLevels = 3
#triangleOptions += "A"
triangleOptions = "q30Dena0.005A"
nLevels = 1


massLumping=False

subgridError = None
shockCapturing = None
subgridError = FFDarcyFC_ASGS(coefficients,nd,stabFlag='2',lag=True)
shockCapturing = ResGradFFDarcy_SC(coefficients,nd,shockCapturingFactor=0.2,lag=False)


multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 10#025
maxLineSearches = 10#0

fullNewtonFlag = True

tolFac = 1.0e-3

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU#PETSc

levelLinearSolver = LU#PETSC
#pick number of layers to use in overlap 
#nLayersOfOverlapForParallel = 1
#parallelPartitioningType = MeshParallelPartitioningTypes.node
#parallelPartitioningType = MeshParallelPartitioningTypes.element
#numericalFluxType = DarcyFCFF_IIPG_exterior

linTolFac = 0.0001

conservativeFlux = {0:'pwl',1:'pwl'}
