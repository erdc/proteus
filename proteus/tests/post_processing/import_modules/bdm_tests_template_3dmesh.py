from proteus import iproteus as ip
from proteus import default_p as p
from proteus import default_n as n
from proteus import default_s,default_so
import numpy
import proteus as pr
from importlib import reload
reload(p)
reload(n)

p.nd = 3
p.name = "BDM2_Test_File"

p.x0 = (0.,0.,0.)
p.L=(1.0,1.0,1.0);
p.T=1.0

p.nc = 1

class velEx(object):
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A = numpy.reshape(self.aex(X),(3,3))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)

def A(x):
    return numpy.array([[1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]],'d')
def f(x):
    return 1.0

class uEx(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2 + x[1]**2 + x[2]**2
    def uOfXT(self,x,T):
        return self.uOfX(x)
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:3],(3,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

eps=1.0e-4
    
def onDirichletBoundary(x):
    if (x[0] <= p.x0[0] + eps or
        x[1] <= p.x0[1] + eps or
        x[1] >= p.x0[1] + p.L[1] - eps or
        x[2] <= p.x0[2] + eps or
        x[2] >= p.x0[2] + p.L[2] - eps):
        return True
    else:
        return False
    
def getDBC(x,flag):
    if onDirichletBoundary(x):
        return lambda x,t: uEx().uOfXT(x,t)

def getAdvFluxBC(x,flag):
    pass

def getDiffFluxBC(x,flag):
    if x[0] >= p.x0[0] + p.L[0] - eps:
        n = numpy.zeros((p.nd,),'d'); n[0]=1.0
        return lambda x,t: numpy.dot(velEx(uEx(),A).uOfXT(x,t),n)
    elif not onDirichletBoundary(x):
        return lambda x,t: 0.0

p.analyticalSolution = {0:uEx()}
p.dirichletConditions = {0:getDBC}
aOfX = {0:A}; fOfX = {0:f}
p.advectiveFluxBoundaryConditions = {0:getAdvFluxBC}
p.diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC}}
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
n.nn = 2
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

n.cfluxtag = 'pwl-bdm2'
n.conservativeFlux = {0:'pwl-bdm2'}

#########################################################################

so = default_so
so.name = p.name 
so.sList=[default_s]

########################################################################
from proteus import *
opts = None
ns = NumericalSolution.NS_base(so,[p],[n],so.sList,ip.opts)
