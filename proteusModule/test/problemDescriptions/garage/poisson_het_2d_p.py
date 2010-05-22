from pyadh import *
from pyadh.default_p import *
"""
Heterogenous Poisson's equation in 2D
"""

##\page Tests Test Problems 
# \ref poisson_het_2d_p.py "Heterogeneous Poisson's equation"
#

##\ingroup test
#\file poisson_het_2d_p.py
#
#\brief Heterogeneous Poisson's equation in 2D

nd = 2

initialConditions = None

Ident = Numeric.zeros((nd,nd),Numeric.Float)
Ident[0,0]=1.0; Ident[1,1] = 1.0

class velEx:
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = Numeric.reshape(self.aex(X),(2,2))
        #mwf debug
        #print """du =%s a= %s v= %s """ % (du,A,-Numeric.dot(A,du))
        return -Numeric.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)


##################################################
#u = 1-x_{i}, a = I, f = 0 
#dir bc's along x, no flow along other axis
ax = 0 #0 1 
def a0(x):
    return Ident.flat[:]
def f0(x):
    return 0.0

class u0Ex:
    def __init__(self):
        self.ax = ax
    def uOfX(self,X):
        return 1.0-X[self.ax]
    def uOfXT(self,X,T):
        return 1.0-X[self.ax]
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.zeros((2,),Numeric.Float)
        du[ax] = -1.0
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)
    
def getDBC0(x):
    if x[ax] == 0.0:
        return lambda x,t: 1.0
    if x[ax] == 1.0:
        return lambda x,t: 0.0

##################################################
##################################################
#u = \sum x_{i}, a = I, f = 0 
#dirichlet bc's everywhere
def a1(x):
    return Ident.flat[:]
def f1(x):
    return 0.0

class u1Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return X[0]+X[1]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.ones((2,),Numeric.Float)
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC1(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u1Ex().uOfXT(x,t)
##################################################
##################################################
#u = \sum x_{i}^2, a = I, f = -2*d 
#dirichlet bc's everywhere
def a2(x):
    return Ident.flat[:]
def f2(x):
    return -2.0*nd

class u2Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return X[0]**2+X[1]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*Numeric.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC2(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u2Ex().uOfXT(x,t)
##################################################
##################################################
#u = sin^2(2\pi x) + cos^2(2\pi y) + x+y + 5.0, a = I,
#f = -8\pi^2cos^2(2\pi x) + 8\pi^2cos^2(2\pi y)
#    +8\pi^2sin^2(2\pi x) - 8\pi^2sin^2(2\pi y)
#dirichlet bc's everywhere

def a3(x):
    return Ident.flat[:]
def f3(x):
    return 8.0*pi**2*(-cos(2.*pi*x[0])**2 + cos(2.0*pi*x[1])**2 + sin(2.0*pi*x[0])**2 - sin(2.0*pi*x[1])**2)

class u3Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return sin(2.0*pi*X[0])**2 + cos(2.0*pi*X[1])**2 + X[0]+X[1]+5.0
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = [1.0+4.0*pi*cos(2.0*pi*X[0])*sin(2.0*pi*X[0]),
              1.0-4.0*pi*cos(2.0*pi*X[1])*sin(2.0*pi*X[1])]
        return Numeric.reshape(Numeric.array(du),(2,))
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC3(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u3Ex().uOfXT(x,t)
##################################################
##################################################
#u = x + y, a_00 = sin^2(2\pi x) + 5, a_11 = cos^2(2\pi y) + 5.0
#f = 4\pi*(-cos(2\pi*x)sin(2\pi*x) + cos(2\pi*y)sin(2\pi*y))
#dirichlet bc's everywhere

def a4(x):
    return Numeric.array([[sin(2.0*pi*x[0])**2 + 5.0,0.0],[0.0,cos(2.0*pi*x[1])**2 + 5.0]],Numeric.Float)
def f4(x):
    return 4.0*pi*(-cos(2.0*pi*x[0])*sin(2.0*pi*x[0]) + cos(2.0*pi*x[1])*sin(2.0*pi*x[1]))

class u4Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]+x[1]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.ones((2,),Numeric.Float)
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC4(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u4Ex().uOfXT(x,t)
##################################################
##################################################
#u = x*x + y*y, a_00 = x + 5, a_11 = y + 5.0
#f = -2*x -2*(5+x) -2*y-2*(5+y)
#dirichlet bc's everywhere

def a5(x):
    return Numeric.array([[x[0] + 5.0,0.0],[0.0,x[1] + 5.0]],Numeric.Float)
def f5(x):
    return -2.0*x[0] -2*(5.+x[0]) -2.*x[1]-2.*(5.+x[1])

class u5Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*Numeric.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC5(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u5Ex().uOfXT(x,t)
##################################################
##################################################
#u = sin^2(2\pi x) + cos^2(2\pi y) + x+y + 5.0,
#a_00 = x*x + 5.0 , a_11 = y*y + 5.0
#f = 2.*x*(1. + 4*pi*Cos(2*pi*x)*Sin(2*pi*x)) - 
#   (5 + x*x)*(8*pi*pi*Cos^2(2*pi*x) - 8*pi*pi*Sin^2(2*pi*x)) - 
#   2.*y*(1. - 4*pi*Cos(2*pi*y)*Sin(2*pi*y)) - 
#   (5 + 1.*y*y)*(-8*pi*pi*Cos^2(2*pi*y) + 8*pi*pi*Sin^2(2*pi*y))

    
#dirichlet bc's everywhere

def a6(x):
    return Numeric.array([[x[0]**2 + 5.0,0.0],[0.0,x[1]**2 + 5.0]],Numeric.Float)
def f6(x):
    return -2.*x[0]*(1. + 4.*pi*cos(2.*pi*x[0])*sin(2.*pi*x[0])) - \
           (5.0 + x[0]**2)*8.*pi*pi*(cos(2.0*pi*x[0])**2 - sin(2.0*pi*x[0])**2) - \
           2.*x[1]*(1. - 4.*pi*cos(2.*pi*x[1])*sin(2.*pi*x[1]))  - \
           (5. + x[1]**2)*8.*pi*pi*(-cos(2.0*pi*x[1])**2 + sin(2.*pi*x[1])**2)

class u6Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return sin(2.0*pi*X[0])**2 + cos(2.0*pi*X[1])**2 + X[0]+X[1]+5.0
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = [1.0+4.0*pi*cos(2.0*pi*X[0])*sin(2.0*pi*X[0]),
              1.0-4.0*pi*cos(2.0*pi*X[1])*sin(2.0*pi*X[1])]
        return Numeric.reshape(Numeric.array(du),(2,))
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC6(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u6Ex().uOfXT(x,t)
##################################################


class PoissonEquationCoefficients(TransportCoefficients.TC_base):
    """
    Variable coefficient Poisson's equation.
    """
    def __init__(self,aOfX,fOfX,nc=1):
        self.aOfX = aOfX
        self.fOfX = fOfX
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        for i in range(nc):
            diffusion[i] = {i : {i:'constant'}}
            reaction[i]  = {i : 'constant'}
            potential[i] = {i : 'u'}
        #end i
        TransportCoefficients.TC_base.__init__(self,
                                               nc,
                                               mass,
                                               advection,
                                               diffusion,
                                               potential,
                                               reaction,
                                               hamiltonian)
    def evaluate(self,t,c):
        for ci in range(self.nc):
            for i in range(len(c[('r',ci)].ravel())):
                c[('r',ci)].flat[i] = -self.fOfX[ci](c['x'].flat[3*i:3*(i+1)])
                c[('a',ci,ci)].flat[nd*nd*i:nd*nd*(i+1)] = self.aOfX[ci](c['x'].flat[3*i:3*(i+1)]).ravel()
            #end i
        #end ci
    #end def


nc = 1

#ex 0
#analyticalSolution = {0:u0Ex(),1:u0Ex()}
#dirichletConditions = {0:getDBC0,1:getDBC0}
#aOfX = {0:a0,1:a0}; fOfX = {0:f0,1:f0}
#ex 1
#analyticalSolution = {0:u1Ex(),1:u1Ex()}
#dirichletConditions = {0:getDBC1,1:getDBC1}
#aOfX = {0:a1,1:a1}; fOfX = {0:f1,1:f1}
#ex 2
#analyticalSolution = {0:u2Ex(),1:u2Ex()}
#dirichletConditions = {0:getDBC2,1:getDBC2}
#aOfX = {0:a2,1:a2}; fOfX = {0:f2,1:f2}
#ex 3
#analyticalSolution = {0:u3Ex(),1:u3Ex()}
#dirichletConditions = {0:getDBC3,1:getDBC3}
#aOfX = {0:a3,1:a3}; fOfX = {0:f3,1:f3}
#ex 4
#analyticalSolution = {0:u4Ex(),1:u4Ex()}
#dirichletConditions = {0:getDBC4,1:getDBC4}
#aOfX = {0:a4,1:a4}; fOfX = {0:f4,1:f4}
#ex 5
#analyticalSolution = {0:u5Ex(),1:u5Ex()}
#dirichletConditions = {0:getDBC5,1:getDBC5}
#aOfX = {0:a5,1:a5}; fOfX = {0:f5,1:f5}
#ex 6
analyticalSolution = {0:u6Ex(),1:u6Ex()}
dirichletConditions = {0:getDBC6,1:getDBC6}
aOfX = {0:a6,1:a6}; fOfX = {0:f6,1:f6}

analyticalSolutionVelocity = {0:velEx(analyticalSolution[0],aOfX[0])}

coefficients = PoissonEquationCoefficients(aOfX,fOfX,nc)
   

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

