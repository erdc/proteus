from pyadh import *
from pyadh.default_n import *
from redist_ellY_2d_p import *

timeIntegration = NoIntegration
nDTout = 1
DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=3
#nLevels = 6
nLevels = 3


multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver

matrix = SparseMatrix

