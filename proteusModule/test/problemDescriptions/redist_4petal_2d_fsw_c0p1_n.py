from pyadh import *
from pyadh.default_n import *
from redist_4petal_2d_p import *

timeIntegration = NoIntegration

DT = None
nDTout = 1


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=5
nLevels = 5

multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FSWEikonalSolver

nl_atol_res = 1.0e-12
tolFac = 1.0e-12
maxNonlinearIts = 100

matrix = SparseMatrix

