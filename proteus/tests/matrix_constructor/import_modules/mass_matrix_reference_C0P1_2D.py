from proteus import iproteus as ip
from proteus import default_p as p
from proteus import default_n as n
from proteus import default_s,default_so
import numpy
from importlib import reload
reload(p)
reload(n)

import proteus as pr

p.nd = 2
p.name = "Mass_matrix_test"

p.rdomain = pr.Domain.unitSimplex(2)
p.polyfile = "reference_triangle_2d"
p.rdomain.writePoly(p.polyfile)
p.domain = None
n.triangleOptions = "Yp"

p.nc = 3

def getDBC(x,flag):
    if x[0] in [0.0] and x[1] in [0.0]:
        pass

p.dirichletConditions = {0:getDBC, 1:getDBC, 2:getDBC}
p.advectiveFluxBoundaryConditions = {}
p.diffusiveFluxBoundaryConditions = {}
p.periodicDirichletConditions = None

p.coefficients = pr.TransportCoefficients.DiscreteMassMatrix(p.nd)

############################

n.timeIntegration = pr.TimeIntegration.NoIntegration
n.nDTout = 1
n.T = 1
n.parallel = False

n.femSpaces = dict((i,pr.FemTools.C0_AffineLinearOnSimplexWithNodalBasis) for i in range(p.nc))
n.elementQuadrature = pr.Quadrature.SimplexGaussQuadrature(p.nd,4)
n.elementBoundaryQuadrature = pr.Quadrature.SimplexGaussQuadrature(p.nd-1,4)
n.nn = 3
n.nLevels = 1

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
