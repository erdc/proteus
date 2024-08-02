from proteus import *
from proteus.default_p import *
"""
Heterogeneous Poisson's equation, -div(a(x)u) = f(x), on unit domain [0,1]x[0,1]x[0,1]
"""

##\page Tests Test Problems
# \ref poisson_3d_p.py "Heterogeneous Poisson's equation, -div(a(x)u) = f(x), on unit domain [0,1]x[0,1]x[0,1]"
#

##\ingroup test
#\file poisson_3d_p.py
#
#\brief Heterogenous Poisson's equations in 3D unit domain [0,1]x[0,1]x[0,1]

#space dimension
nd = 3
#if unstructured would need variable polyfile or meshfile set
x0 = (-3.,-3.,-3.)
L  = ( 6., 6., 6.)
domain = None
polyfile = None

test_hexMesh_3x3 = False
if test_hexMesh_3x3 == True:
    meshfile='hexMesh_3x3'
    domain = Domain.MeshHexDomain(meshfile)
    x0 = (-3.,-3.,-3.)
    L  = ( 6., 6., 6.)
#steady-state so no initial conditions
initialConditions = None
#use sparse diffusion representation
sd=True
#identity tensor for defining analytical heterogeneity functions
Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0; Ident[2,2]=1.0


#for computing exact 'Darcy' velocity
class velEx(object):
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = numpy.reshape(self.aex(X),(3,3))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)


##################################################
#define coefficients a(x)=[a_{ij}] i,j=0,2, right hand side f(x)  and analytical solution u(x)
#u = x*x + y*y + z*z, a_00 = x + 5, a_11 = y + 5.0 + a_22 = z + 10.0
#f = -2*x -2*(5+x) -2*y-2*(5+y) -2*z-2*(10+z)
#

def a5(x):
    return numpy.array([[x[0] + 5.0,0.0,0.0],[0.0,x[1] + 5.0,0.0],[0.0,0.0,x[2]+10.0]],'d')
def f5(x):
    return -2.0*x[0] -2*(5.+x[0]) -2.*x[1]-2.*(5.+x[1]) -2.*x[2]-2.*(10+x[2])
#'manufactured' analytical solution
class u5Ex(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2+x[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:3],(3,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

eps=1.0e-4

def onDirichletBoundary(x):
    if (x[0] <= x0[0] + eps or
        x[1] <= x0[1] + eps or
        x[1] >= x0[1] + L[1] - eps or
        x[2] <= x0[2] + eps or
        x[2] >= x0[2] + L[2] - eps):
        return True
    else:
        return False

#dirichlet boundary condition functions on (x=0,y,z), (x,y=0,z), (x,y=1,z), (x,y,z=0), (x,y,z=1)
def getDBC5(x,flag):
    if onDirichletBoundary(x):
        return lambda x,t: u5Ex().uOfXT(x,t)
def getAdvFluxBC5(x,flag):
    return None

#specify flux on (x=1,y,z)
def getDiffFluxBC5(x,flag):
    if x[0] >= x0[0] + L[0] - eps:
        n = numpy.zeros((nd,),'d');
        n[0]=1.0
        return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
    elif not onDirichletBoundary(x):
        return lambda x,t: 0.0
#store a,f in dictionaries since coefficients class allows for one entry per component
aOfX = {0:a5}; fOfX = {0:f5}

#one component
nc = 1
#load analytical solution, dirichlet conditions, flux boundary conditions into the expected variables
analyticalSolution = {0:u5Ex()}
analyticalSolutionVelocity = {0:velEx(analyticalSolution[0],aOfX[0])}
#
dirichletConditions = {0:getDBC5}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC5}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC5}}

#equation coefficient names
coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,nc,nd)
#
coefficients.variableNames=['u0']
