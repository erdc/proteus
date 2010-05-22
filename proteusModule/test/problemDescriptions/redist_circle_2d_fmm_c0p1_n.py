from pyadh import *
from pyadh.default_n import *
from redist_circle_2d_p import *

timeIntegration = NoIntegration
nDTout = 1
DT = None

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn=11#3
#nLevels = 6
nLevels = 3


multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver

matrix = SparseMatrix

