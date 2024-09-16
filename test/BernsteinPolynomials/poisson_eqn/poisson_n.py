from proteus import *
from proteus.default_n import *
try:
    from .parameters_for_poisson import *
    from .poisson_p import *
except:
    from parameters_for_poisson import *
    from poisson_p import *

parallel = False
#0 there is no special treatment for time derivative
timeIntegration = NoIntegration
quad_order=2*ct.pDegree+1
if ct.useHex:
    domain.MeshOptions.hex=True
    domain.MeshOptions.quad=True
    if ct.pDegree==1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif ct.pDegree==2:
        if ct.useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % ct.pDegree_vof)
    elementQuadrature = CubeGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,quad_order)
else:
    domain.MeshOptions.hex=False
    domain.MeshOptions.quad=False
    hex = False
    quad = False
    if ct.pDegree == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif ct.pDegree == 2:
        if ct.useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnSimplex}
        else:
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % ct.pDegree_vof)
    elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)


nLevels = 1 
numericalFluxType=None
subgridError = None
shockCapturing = None

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = Newton
#2 since problem is linear
maxNonlinearIts = 1

fullNewtonFlag = False
tolFac = 0.0
nl_atol_res = 1.0e-12

matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
linearSolverConvergenceTest= 'r'

linTolFac = 0.0 
l_atol_res = 1.0e-12
