from pyadh import *
from pyadh.default_n import *
from redist_4petal_2d_p import *

timeIntegration = NoIntegration

DT = None
nDTout = 1


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=51
nLevels = 1#5

multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver


matrix = SparseMatrix

