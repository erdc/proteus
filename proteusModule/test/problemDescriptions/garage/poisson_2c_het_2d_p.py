from pyadh import *
from pyadh.default_p import *
"""
Heterogeneous Poisson's equations for two components (uncoupled) in 2D
"""

##\page Tests Test Problems 
# \ref poisson_2c_het_2d_p.py "Heterogeneous Poisson's equation for two components"
#

##\ingroup test
#\file poisson_2c_het_2d_p.py
#
#\brief Heterogensou Poisson's equations for two components (uncoupled) in 2D

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

def getAdvFluxBC0(x):
    pass
def getDiffFluxBC0(x):
    pass

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
    if x[0] in [1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u1Ex().uOfXT(x,t)

def getAdvFluxBC1(x):
    pass
def getDiffFluxBC1(x):
    if x[0] == 0.0:
        n = Numeric.array([-1.0,0.0],Numeric.Float)
        #mwf debug
        #print """getDiffFlux1 x= %s q.n= %s """ % (x,Numeric.dot(velEx(u1Ex(),a1).uOfXT(x,0.0),n))
        return lambda x,t: Numeric.dot(velEx(u1Ex(),a1).uOfXT(x,t),n)

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
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u2Ex().uOfXT(x,t)
def getAdvFluxBC2(x):
    pass
def getDiffFluxBC2(x):
    if x[0] == 1.0:
        n = Numeric.zeros((nd,),Numeric.Float); n[0]=1.0
        #mwf debug
        #print """getDiffFlux2 x= %s q.n= %s """ % (x,Numeric.dot(velEx(u2Ex(),a2).uOfXT(x,0.0),n))        
        return lambda x,t: Numeric.dot(velEx(u2Ex(),a2).uOfXT(x,t),n)


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
    if x[0] in [0.0,1.0] or x[1] in [0.0]:
        return lambda x,t: u3Ex().uOfXT(x,t)
def getAdvFluxBC3(x):
    pass
def getDiffFluxBC3(x):
    if x[1] == 1.0:
        n = Numeric.zeros((nd,),Numeric.Float); n[1]=1.0
        return lambda x,t: Numeric.dot(velEx(u3Ex(),a3).uOfXT(x,t),n)
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
    if x[0] in [0.0,1.0] or x[1] in [1.0]:
        return lambda x,t: u4Ex().uOfXT(x,t)
#for pwl might be easier to enforce boundary flux as advective flux
def getAdvFluxBC4(x):
    if x[1] == 0.0:
        n = Numeric.zeros((nd,),Numeric.Float); n[1]=-1.0
        return lambda x,t: Numeric.dot(velEx(u4Ex(),a4).uOfXT(x,t),n)
def getDiffFluxBC4(x):
    pass
#def getAdvFluxBC4(x):
#    pass
#def getDiffFluxBC4(x):
#     if x[1] == 0.0:
#         n = Numeric.zeros((nd,),Numeric.Float); n[1]=-1.0
#         return lambda x,t: Numeric.dot(velEx(u4Ex(),a4).uOfXT(x,t),n)
    
##################################################
##################################################
#u = x*x + y*y, a_00 = x + 5, a_11 = y + 5.0
#f = -2*x -2*(5+x) -2*y-2*(5+y)
#

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
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u5Ex().uOfXT(x,t)
def getAdvFluxBC5(x):
    if x[0] == 1.0:
        n = Numeric.zeros((nd,),Numeric.Float); n[0]=1.0
        return lambda x,t: Numeric.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
def getDiffFluxBC5(x):
    pass
#def getAdvFluxBC5(x):
#    pass
#def getDiffFluxBC5(x):
#    if x[0] == 1.0:
#        n = Numeric.zeros((nd,),Numeric.Float); n[0]=1.0
#        return lambda x,t: Numeric.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
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
eps = 1.0e-8
def getDBC6(x):
    #if x[0] in [0.0,1.0] or x[1] in [0.0]:
    if abs(x[0]-0.0) < eps or abs(x[0]-1.0) < eps or abs(x[1]-0.0) < eps:
        return lambda x,t: u6Ex().uOfXT(x,t)
def getAdvFluxBC6(x):
    if abs(x[1]-1.0) < eps:
        n = Numeric.zeros((nd,),Numeric.Float); n[1]=1.0
        return lambda x,t: Numeric.dot(velEx(u6Ex(),a6).uOfXT(x,t),n)
def getDiffFluxBC6(x):
    pass
# def getAdvFluxBC6(x):
#     pass
# def getDiffFluxBC6(x):
#     if x[1] == 1.0:
#         n = Numeric.zeros((nd,),Numeric.Float); n[1]=1.0
#         return lambda x,t: Numeric.dot(velEx(u6Ex(),a6).uOfXT(x,t),n)
##################################################
##################################################
#u = x*x + y*y, a_00 = x*x + 5, a_11 = y*y + 5.0
#f = -(6*x*x + 6*y*y + 20)
#

def a7(x):
    return Numeric.array([[x[0]**2 + 5.0,0.0],[0.0,x[1]**2 + 5.0]],Numeric.Float)
def f7(x):
    return -6.0*x[0]**2 -6.0*x[1]**2 -20.0

class u7Ex:
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

def getDBC7(x):
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u7Ex().uOfXT(x,t)
def getAdvFluxBC7(x):
    if x[0] == 1.0:
        n = Numeric.zeros((nd,),Numeric.Float); n[0]=1.0
        return lambda x,t: Numeric.dot(velEx(u7Ex(),a7).uOfXT(x,t),n)
def getDiffFluxBC7(x):
    pass
##################################################
##################################################
#u = x*x + y*y, a_00 = x*x*y + 5, a_11 = y*y*x + 5.0
#f = -(6*x*x*y + 6*y*y*x + 20)
#

def a8(x):
    return Numeric.array([[x[1]*x[0]**2 + 5.0,0.0],[0.0,x[0]*x[1]**2 + 5.0]],Numeric.Float)
def f8(x):
    return -6.0*x[1]*x[0]**2 -6.0*x[0]*x[1]**2 -20.0

class u8Ex:
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

def getDBC8(x):
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u8Ex().uOfXT(x,t)
def getAdvFluxBC8(x):
    if x[0] == 1.0:
        n = Numeric.zeros((nd,),Numeric.Float); n[0]=1.0
        return lambda x,t: Numeric.dot(velEx(u8Ex(),a8).uOfXT(x,t),n)
def getDiffFluxBC8(x):
    pass

# def getDBC8(x):
#     if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
#         return lambda x,t: u8Ex().uOfXT(x,t)
# def getAdvFluxBC8(x):
#     pass
# def getDiffFluxBC8(x):
#     pass
##################################################


nc = 2

#ex 0
#analyticalSolution = {0:u0Ex(),1:u0Ex()}
#dirichletConditions = {0:getDBC0,1:getDBC0}
#aOfX = {0:a0,1:a0}; fOfX = {0:f0,1:f0}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC0,1:getAdvFluxBC0}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC0},1:{1:getDiffFluxBC0}}
#
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC0,1:getAdvFluxBC0}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC0},1:{1:getDiffFluxBC0}}
#
#ex 1
#analyticalSolution = {0:u1Ex(),1:u1Ex()}
#dirichletConditions = {0:getDBC1,1:getDBC1}
#aOfX = {0:a1,1:a1}; fOfX = {0:f1,1:f1}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC1,1:getAdvFluxBC1}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC1},1:{1:getDiffFluxBC1}}
#
#ex 2
#analyticalSolution = {0:u2Ex(),1:u2Ex()}
#dirichletConditions = {0:getDBC2,1:getDBC2}
#aOfX = {0:a2,1:a2}; fOfX = {0:f2,1:f2}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC2,1:getAdvFluxBC2}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC2},1:{1:getDiffFluxBC2}}
#
#ex 3
#analyticalSolution = {0:u3Ex(),1:u3Ex()}
#dirichletConditions = {0:getDBC3,1:getDBC3}
#aOfX = {0:a3,1:a3}; fOfX = {0:f3,1:f3}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC3,1:getAdvFluxBC3}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC3},1:{1:getDiffFluxBC3}}
#
#ex 4
#analyticalSolution = {0:u4Ex(),1:u4Ex()}
#dirichletConditions = {0:getDBC4,1:getDBC4}
#aOfX = {0:a4,1:a4}; fOfX = {0:f4,1:f4}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC4,1:getAdvFluxBC4}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC4},1:{1:getDiffFluxBC4}}
#
#ex 5
#analyticalSolution = {0:u5Ex(),1:u5Ex()}
#dirichletConditions = {0:getDBC5,1:getDBC5}
#aOfX = {0:a5,1:a5}; fOfX = {0:f5,1:f5}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC5,1:getAdvFluxBC5}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC5},1:{1:getDiffFluxBC5}}
#
#ex 6
analyticalSolution = {0:u6Ex(),1:u6Ex()}
dirichletConditions = {0:getDBC6,1:getDBC6}
aOfX = {0:a6,1:a6}; fOfX = {0:f6,1:f6}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC6,1:getAdvFluxBC6}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC6},1:{1:getDiffFluxBC6}}

#mix ex0 and ex6 
#analyticalSolution = {0:u0Ex(),1:u6Ex()}
#dirichletConditions = {0:getDBC0,1:getDBC6}
#aOfX = {0:a0,1:a6}; fOfX = {0:f0,1:f6}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC0,1:getAdvFluxBC6}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC0},1:{1:getDiffFluxBC6}}

#ex 7
#analyticalSolution = {0:u7Ex(),1:u7Ex()}
#dirichletConditions = {0:getDBC7,1:getDBC7}
#aOfX = {0:a7,1:a7}; fOfX = {0:f7,1:f7}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC7,1:getAdvFluxBC7}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC7},1:{1:getDiffFluxBC7}}

#ex8
#analyticalSolution = {0:u8Ex(),1:u8Ex()}
#dirichletConditions = {0:getDBC8,1:getDBC8}
#aOfX = {0:a8,1:a8}; fOfX = {0:f8,1:f8}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC8,1:getAdvFluxBC8}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC8},1:{1:getDiffFluxBC8}}

analyticalSolutionVelocity = {0:velEx(analyticalSolution[0],aOfX[0]),
                              1:velEx(analyticalSolution[1],aOfX[1])}

coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,nc,nd)

coefficients.variableNames=['u0','u1']

   
fluxBoundaryConditions = {0:'setFlow',1:'setFlow'}



