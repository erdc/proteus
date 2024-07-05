from proteus.iproteus import *
from importlib import reload
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
p.name = "TwoPhase_mass_matrix_test"

p.x0 = [-1.0,-1.0]
p.L = [2,2]
p.nnx=p.nny=3

p.nc = 3

def getDBC(x,flag):
    if x[0] in [0.0] and x[1] in [0.0]:
        pass

p.dirichletConditions = {0:getDBC, 1:getDBC, 2:getDBC}
p.advectiveFluxBoundaryConditions = {}
p.diffusiveFluxBoundaryConditions = {}
p.periodicDirichletConditions = None

phase_func = lambda x: x[0]
p.coefficients = pr.TransportCoefficients.DiscreteMassMatrix(nd = p.nd)

############################

n.timeIntegration = pr.TimeIntegration.NoIntegration
n.nDTout = 1
n.T = 1
n.parallel = False

n.femSpaces = {0:pr.FemTools.C0_AffineLinearOnCubeWithNodalBasis,
               1:pr.FemTools.Q2,
               2:pr.FemTools.Q2}

n.quad = True
n.elementQuadrature = pr.Quadrature.CubeGaussQuadrature(p.nd,4)
n.elementBoundaryQuadrature = pr.Quadrature.CubeGaussQuadrature(p.nd-1,4)
n.nn = 3
n.nLevels = 1

#########################################################################

so = default_so
so.name = p.name 
so.sList=[default_s]

########################################################################
from proteus import *
opts = None
ns = NumericalSolution.NS_base(so,[p],[n],so.sList,ip.opts)