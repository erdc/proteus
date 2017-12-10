from proteus import *
from proteus.default_n import *
from sm_gl_6_3d_p import *

timeIntegration = NoIntegration
stepControl = Newton_controller
#tnList=[0,T,2*T]#,3*T]
tnList=[float(i) for i in range(2)]
#tnList=[0,T]

elementQuadrature = SimplexGaussQuadrature(nd,6)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,6)

#elementQuadrature = SimplexGaussQuadrature(nd,6)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,6)

femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
             1:C0_AffineQuadraticOnSimplexWithNodalBasis,
             2:C0_AffineQuadraticOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
#femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#             1:C0_AffineLinearOnSimplexWithNodalBasis,
#             2:C0_AffineLinearOnSimplexWithNodalBasis}
subgridError = None

massLumping = False

numericalFluxType = Stress_IIPG_exterior#None
#numericalFluxType = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True
maxNonlinearIts=25
maxLineSearches =25
maxSolverFailures = 1

tolFac = 0.0

atol = 1.0e-4
nl_rtol_res =0.0
nl_atol_res =1.0e-4
matrix = SparseMatrix

#multilevelLinearSolver = PETSc#LU
multilevelLinearSolver = LU

#levelLinearSolver = PETSc#LU
levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.0

conservativeFlux = None

#auxiliaryVariables=[PlasticWork()]
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1
#parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
#nLayersOfOverlapForParallel = 0
