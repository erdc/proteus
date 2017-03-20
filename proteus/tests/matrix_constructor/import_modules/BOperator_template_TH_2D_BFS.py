from proteus.iproteus import *
reload(default_p)
reload(default_n)
reload(default_so)
from proteus import iproteus as ip
from proteus import default_p as p
from proteus import default_n as n
from proteus import default_s,default_so
import numpy
import proteus as pr

p.nd = 2
p.name = "Laplace_matrix_test"

p.x0 = [-1.0,-1.0]
p.L = [2,2]

p.domain = Domain.PlanarStraightLineGraphDomain(
    vertices   =[[-1.0,0.0],[-1.0,1.0],[5.0,1.0],[5.0,-1.0],[0.0,-1.0],[0.0,0.0]],
    vertexFlags=[ 1        , 1        , 2       , 4        , 4        , 4       ],
    segments    =[[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]],
    segmentFlags=[ 1   , 2   , 3   , 4   , 4   , 4   ],
    name="backwardsFacingStep",
    units='cm')
p.domain.writePoly('backwardsFacingStep')
p.he = 4.0

p.nc = 3

def getDBC(x,flag):
    if x[0] in [0.0] and x[1] in [0.0]:
        pass

p.dirichletConditions = {0:getDBC, 1:getDBC, 2:getDBC}
p.advectiveFluxBoundaryConditions = {}
p.diffusiveFluxBoundaryConditions = {}
p.periodicDirichletConditions = None

p.coefficients = pr.TransportCoefficients.Stokes(p.nd)

############################

n.timeIntegration = pr.TimeIntegration.NoIntegration
n.nDTout = 1
n.T = 1
n.parallel = False

n.femSpaces = {0:pr.FemTools.C0_AffineLinearOnSimplexWithNodalBasis,
               1:pr.FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis,
               2:pr.FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis}

n.elementQuadrature = pr.Quadrature.SimplexGaussQuadrature(p.nd,4)
n.elementBoundaryQuadrature = pr.Quadrature.SimplexGaussQuadrature(p.nd-1,4)
n.nn = 3
n.nLevels = 1

n.subgridError = None
n.shockCapturing = None
n.multilevelNonlinearSolver = pr.NonlinearSolvers.Newton
n.levelNonlinearSolver = pr.NonlinearSolvers.Newton
n.maxNonlinearIts = 1
n.fullNewtonFlag = True
n.totFac = 1.0e-8
n.nl_atol_res = 1.0e-8
n.matrix = pr.LinearAlgebraTools.SparseMatrix
n.multilevelLinearSolver = pr.LinearSolvers.LU
n.levelLinearSolver = pr.LinearSolvers.LU#MGM#
n.linearSolverConvergenceTest= 'r'#r-true'#'r'

#########################################################################

so = default_so
so.name = p.name 
so.sList=[default_s]

########################################################################
from proteus import *
opts = None
ns = NumericalSolution.NS_base(so,[p],[n],so.sList,ip.opts)
