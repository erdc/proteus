from builtins import range
from builtins import object
from proteus import iproteus as ip
from proteus import default_p as p
from proteus import default_n as n
from proteus import default_s,default_so
import numpy
import proteus as pr

p.nd = 2
p.name = "BDM2_Test_File"

p.L=(1.0,1.0); 
p.T=1.0

p.nc = 1

class velEx(object):
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A = numpy.reshape(self.aex(X),(2,2))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)

def A(x):
    return numpy.array([[1.0, 0.0],[0.0, 1.0]],'d')
def f(x):
    return 0.0

class uEx(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]*x[1]
    def uOfXT(self,x,T):
        return self.uOfX(x)
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: uEx().uOfXT(x,t)

def getAdvFluxBC(x,flag):
    pass

def getDiffFluxBC(x,flag):
    pass
 
p.analyticalSolution = {0:uEx()}
p.dirichletConditions = {0:getDBC}
aOfX = {0:A}; fOfX = {0:f}
p.advectiveFluxBoundaryConditions = {0:getAdvFluxBC}
#p.diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC}}
p.periodicDirichletConditions = None

p.coefficients = pr.TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,p.nc,p.nd)

############################

n.timeIntegration = pr.TimeIntegration.NoIntegration
n.nDTout = 1
n.T = 1
n.parallel = False

n.femSpaces = dict((i,pr.FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis) for i in range(p.nc))
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

if n.parallel:
    n.multilevelLinearSolver = pr.KSP_petsc4py#PETSc#LU
    n.levelLinearSolver = pr.KSP_petsc4py#PETSc#LU#MGM#PETSc#
    n.nLayersOfOverlapForParallel = 1
    n.parallelPartitioningType = pr.MeshParallelPartitioningTypes.element
    n.numericalFluxType = pr.Advection_DiagonalUpwind_Diffusion_IIPG_exterior
    n.linearSmoother = None
else:
    n.multilevelLinearSolver = pr.LinearSolvers.LU
    n.levelLinearSolver = pr.LinearSolvers.LU#MGM#
    n.linearSolverConvergenceTest= 'r'#r-true'#'r'



n.linearSmoother = pr.LinearSolvers.StarILU#GaussSeidel#Jacobi#StarILU

n.linTolFac = 0.0
n.l_atol_res = 1.0e-10

n.multigridCycles = 0

n.cfluxtag  = 'pwl-bdm2'#'pwl-bdm'#'sun-rt0','sun-gs-rt0','pwc','pwl','pwl-bdm','point-eval'
n.conservativeFlux =  {0:'pwl-bdm2'}#,1:'pwl-ib-fix-0'}#dict((i,cfluxtag) for i in range(nc))
#need this for sun-wheeler-gs
#if n.cfluxtag == 'sun-gs-rt0':
#n.numericalFluxType = pr.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior
#n.useStrongDirichletConstraints = True
#n.preSmooths = 3
#n.postSmooths = 3

#########################################################################

so = default_so
so.name = p.name 
so.sList=[default_s]

########################################################################
from proteus import *
opts = None
simFlagsList = [{}]
simFlagsList[0]['storeQuantities'] = ['u',"q:('velocity',0)"]
ns = NumericalSolution.NS_base(so,[p],[n],so.sList,ip.opts,simFlagsList)
import pdb
pdb.set_trace()

failed = ns.calculateSolution('ladr_run1')
assert(not failed)
