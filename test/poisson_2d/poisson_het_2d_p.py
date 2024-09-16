from proteus import *
from proteus.default_p import *
from proteus import defaults
defaults.reset_default_p()
"""
Heterogeneous Poisson's equations for one component (uncoupled) in 2D
"""

##\page Tests Test Problems 
# \ref poisson_het_2d_p.py "Heterogeneous Poisson's equation for 1 components"
#

##\ingroup test
#\file poisson_het_2d_p.py
#
#\brief Heterogensou Poisson's equations for one component (uncoupled) in 2D

nd = 2

test_problem_number = 5
assert test_problem_number in [2,5]

initialConditions = None
sd=True
Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0
domain = None
polyfile = None

class velEx(object):
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = numpy.reshape(self.aex(X),(2,2))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)


##################################################
##################################################
#u = \sum x_{i}^2, a = I, f = -2*d 
#dirichlet bc's everywhere
def a2(x):
    return Ident.flat[:]
def f2(x):
    return -2.0*nd

class u2Ex(object):
    def __init__(self):
        pass
    def uOfX(self,X):
        return X[0]**2+X[1]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC2(x,flag):
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u2Ex().uOfXT(x,t)
def getDiffFluxBC2(x,flag): 
    if x[0] == 1.0:
        n = numpy.zeros((nd,),'d'); n[0]=1.0
        return lambda x,t: numpy.dot(velEx(u2Ex(),a2).uOfXT(x,t),n)
def getAdvFluxBC2(x,flag):
    pass


##################################################
##################################################
#u = x*x + y*y, a_00 = x + 5, a_11 = y + 5.0
#f = -2*x -2*(5+x) -2*y-2*(5+y)
#

def a5(x):
    return numpy.array([[x[0] + 5.0,0.0],[0.0,x[1] + 5.0]],'d')
def f5(x):
    return -2.0*x[0] -2*(5.+x[0]) -2.*x[1]-2.*(5.+x[1])

class u5Ex(object):
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

eps=1.0e-8
def getDBC5(x,flag):
    if x[0] < eps or x[1] < eps or x[1] > 1.0-eps:
        return lambda x,t: u5Ex().uOfXT(x,t)
def getAdvFluxBC5(x,flag):
   pass
def getDiffFluxBC5(x,flag):
    if x[0] > 1.0-eps:
        n = numpy.zeros((nd,),'d'); n[0]=1.0
        return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
    if flag == 0:
        return lambda x,t: 0.0

#def getAdvFluxBC5(x,flag):
#    pass
#def getDiffFluxBC5(x,flag):
#    if x[0] == 1.0:
#        n = numpy.zeros((nd,),'d'); n[0]=1.0
#        return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)

##################################################

nc = 1

if test_problem_number == 2:
    #ex 2
    analyticalSolution = dict((i,u2Ex()) for i in range(nc))   
    dirichletConditions = dict((i,getDBC2) for i in range(nc)) 
    aOfX = dict((i,a2) for i in range(nc)) ; fOfX = dict((i,f2) for i in range(nc)) 
    advectiveFluxBoundaryConditions =  dict((i,getAdvFluxBC2) for i in range(nc))   
    diffusiveFluxBoundaryConditions = dict((i,{i:getDiffFluxBC2}) for i in range(nc))
else:
    #ex 5
    analyticalSolution = dict((i,u5Ex()) for i in range(nc))  
    dirichletConditions = dict((i,getDBC5) for i in range(nc)) 
    aOfX = dict((i,a5) for i in range(nc)) ; fOfX = dict((i,f5) for i in range(nc)) 
    advectiveFluxBoundaryConditions =  dict((i,getAdvFluxBC5) for i in range(nc))   
    diffusiveFluxBoundaryConditions = dict((i,{i:getDiffFluxBC5}) for i in range(nc))
#

analyticalSolutionVelocity = dict((i,velEx(analyticalSolution[i],aOfX[i])) for i in range(nc))

coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,nc,nd)

coefficients.variableNames=list('u%d' % i for i in range(nc))
   
fluxBoundaryConditions = dict((i,'setFlow') for i in range(nc))


