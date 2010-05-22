from pyadh import *
from pyadh.default_n import *
from redist_circle_3d_p import *

timeIntegration = NoIntegration
nDTout = 1
DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=21
nLevels = 1


multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FSWEikonalSolver

nl_atol_res = 1.0e-8
tolFac = 1.0e-8
maxNonlinearIts = 100

matrix = SparseMatrix

