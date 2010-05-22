from pyadh import *
from pyadh.default_p import *
"""
Heterogeneous Poisson's equations for two components (uncoupled) in 3D
"""

##\page Tests Test Problems 
# \ref poisson_2c_het_3d_p.py "Heterogeneous Poisson's equation for two components"
#

##\ingroup test
#\file poisson_2c_het_3d_p.py
#
#\brief Heterogensou Poisson's equations for two components (uncoupled) in 3D

#mwf for now some test mass conservation depends on mass conservation to
#determine which boundaries are Neumann boundaries

nd = 3

initialConditions = None

Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0; Ident[2,2]=1.0

class velEx:
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = numpy.reshape(self.aex(X),(nd,nd))
        #mwf debug
        #print """du =%s a= %s v= %s """ % (du,A,-numpy.dot(A,du))
        return -numpy.dot(A,du)
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
        du = numpy.zeros((nd,),'d')
        du[ax] = -1.0
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)
    
def getDBC0(x,flag):
    if x[ax] == 0.0:
        return lambda x,t: 1.0
    if x[ax] == 1.0:
        return lambda x,t: 0.0

def getDiffFluxBC0(x,flag):
    pass
def getAdvFluxBC0(x,flag):
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
        return X[0]+X[1]+X[2]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = numpy.ones((nd,),'d')
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC1(x,flag):
    if x[0] in [1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u1Ex().uOfXT(x,t)
def getDiffFluxBC1(x,flag):
    pass
def getAdvFluxBC1(x,flag):
    if x[0] == 0.0:
        n = numpy.array([-1.0,0.0,0.0],'d')
        #mwf debug
        #print """getAdvFlux1 x= %s q.n= %s """ % (x,numpy.dot(velEx(u1Ex(),a1).uOfXT(x,0.0),n))
        return lambda x,t: numpy.dot(velEx(u1Ex(),a1).uOfXT(x,t),n)
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
        return X[0]**2+X[1]**2+X[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:nd],(nd,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC2(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0]:
        return lambda x,t: u2Ex().uOfXT(x,t)
def getDiffFluxBC2(x,flag):
    pass
def getAdvFluxBC2(x,flag):
    if x[2] == 1.0:
        n = numpy.zeros((nd,),'d'); n[2]=1.0
        #mwf debug
        #print """getAdvFlux2 x= %s q.n= %s """ % (x,numpy.dot(velEx(u2Ex(),a2).uOfXT(x,0.0),n))        
        return lambda x,t: numpy.dot(velEx(u2Ex(),a2).uOfXT(x,t),n)
##################################################
##################################################
#u = sin^2(2\pi x) + cos^2(2\pi y) + sin^2(2\pi z)+ x+y+z + 5.0, a = I,
#f = -8\pi^2cos^2(2\pi x) + 8\pi^2cos^2(2\pi y) -8\pi^2cos^2(2\pi z)
#    +8\pi^2sin^2(2\pi x) - 8\pi^2sin^2(2\pi y) +8\pi^2sin^2(2\pi z) 
#dirichlet bc's everywhere

def a3(x):
    return Ident.flat[:]
def f3(x):
    return 8.0*pi**2*(-cos(2.*pi*x[0])**2 + cos(2.0*pi*x[1])**2 + sin(2.0*pi*x[0])**2 - sin(2.0*pi*x[1])**2 -cos(2.*pi*x[2])**2 + sin(2.0*pi*x[2])**2)

class u3Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return sin(2.0*pi*X[0])**2 + cos(2.0*pi*X[1])**2 + sin(2.0*pi*X[2])**2 + X[0]+X[1]+X[2]+5.0
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = [1.0+4.0*pi*cos(2.0*pi*X[0])*sin(2.0*pi*X[0]),
              1.0-4.0*pi*cos(2.0*pi*X[1])*sin(2.0*pi*X[1]),
              1.0+4.0*pi*cos(2.0*pi*X[2])*sin(2.0*pi*X[2])]
        return numpy.reshape(numpy.array(du),(nd,))
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC3(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [0.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u3Ex().uOfXT(x,t)
def getDiffFluxBC3(x,flag):
    pass
def getAdvFluxBC3(x,flag):
    if x[1] == 1.0:
        n = numpy.zeros((nd,),'d'); n[1]=1.0
        return lambda x,t: numpy.dot(velEx(u3Ex(),a3).uOfXT(x,t),n)
##################################################
##################################################
#u = x + y + z, a_00 = sin^2(2\pi x) + 5, a_11 = cos^2(2\pi y) + 5.0
# a_22 = sin^2(2\pi z) + 5
#f = 4\pi*(-cos(2\pi*x)sin(2\pi*x) + cos(2\pi*y)sin(2\pi*y)
#          -cos(2\pi*z)sin(2\pi*z))
#dirichlet bc's everywhere

def a4(x):
    return numpy.array([[sin(2.0*pi*x[0])**2 + 5.0,0.0,0.0],[0.0,cos(2.0*pi*x[1])**2 + 5.0,0.0],
                          [0.0,0.0,sin(2.0*pi*x[2])**2 + 5.0]],
                         'd')
def f4(x):
    return 4.0*pi*(-cos(2.0*pi*x[0])*sin(2.0*pi*x[0]) + cos(2.0*pi*x[1])*sin(2.0*pi*x[1])
                   -cos(2.0*pi*x[2])*sin(2.0*pi*x[2]))

class u4Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]+x[1]+x[2]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = numpy.ones((nd,),'d')
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC4(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u4Ex().uOfXT(x,t)
def getDiffFluxBC4(x,flag):
    pass
def getAdvFluxBC4(x,flag):
    if x[1] == 0.0:
        n = numpy.zeros((nd,),'d'); n[1]=-1.0
        return lambda x,t: numpy.dot(velEx(u4Ex(),a4).uOfXT(x,t),n)
##################################################
##################################################
#u = x*x + y*y + z*z, a_00 = x + 5, a_11 = y + 5.0, a_22 = z + 5.0
#f = -2*x -2*(5+x) -2*y-2*(5+y) -2*z -2*(5+z)
#dirichlet bc's everywhere

def a5(x):
    return numpy.array([[x[0] + 5.0,0.0,0.0],[0.0,x[1] + 5.0,0.0],
                          [0.0,0.0,x[2]+5.0]],'d')
def f5(x):
    return -2.0*x[0] -2*(5.+x[0]) -2.*x[1]-2.*(5.+x[1]) -2.0*x[2] -2*(5.+x[2])

class u5Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2+x[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:nd],(nd,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC5(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0]:
        return lambda x,t: u5Ex().uOfXT(x,t)
def getDiffFluxBC5(x,flag):
    pass
def getAdvFluxBC5(x,flag):
    if x[2] == 1.0:
        n = numpy.zeros((nd,),'d'); n[2]=1.0
        return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
##################################################
##################################################
#u = sin^2(2\pi x) + cos^2(2\pi y) + sin^2(2\pi z) +x+y+z + 5.0,
#a_00 = x*x + 5.0 , a_11 = y*y + 5.0, a_22 = z*z + 5.0
#f = -2.*x*(1. + 4*pi*Cos(2*pi*x)*Sin(2*pi*x)) - 
#   (5 + x*x)*(8*pi*pi*Cos^2(2*pi*x) - 8*pi*pi*Sin^2(2*pi*x)) - 
#   2.*y*(1. - 4*pi*Cos(2*pi*y)*Sin(2*pi*y)) - 
#   (5 + 1.*y*y)*(-8*pi*pi*Cos^2(2*pi*y) + 8*pi*pi*Sin^2(2*pi*y)) -
#   2.*z*(1. + 4*pi*Cos(2*pi*z)*Sin(2*pi*z)) - 
#   (5 + z*z)*(8*pi*pi*Cos^2(2*pi*z) - 8*pi*pi*Sin^2(2*pi*z))


def a6(x):
    return numpy.array([[x[0]**2 + 5.0,0.0,0.0],[0.0,x[1]**2 + 5.0,0.0],
                          [0.0,0.0,x[2]**2 + 5.0]],'d')
def f6(x):
    return -2.*x[0]*(1. + 4.*pi*cos(2.*pi*x[0])*sin(2.*pi*x[0])) - \
           (5.0 + x[0]**2)*8.*pi*pi*(cos(2.0*pi*x[0])**2 - sin(2.0*pi*x[0])**2) - \
           2.*x[1]*(1. - 4.*pi*cos(2.*pi*x[1])*sin(2.*pi*x[1]))  - \
           (5. + x[1]**2)*8.*pi*pi*(-cos(2.0*pi*x[1])**2 + sin(2.*pi*x[1])**2) - \
           2.*x[2]*(1. + 4.*pi*cos(2.*pi*x[2])*sin(2.*pi*x[2])) - \
           (5.0 + x[2]**2)*8.*pi*pi*(cos(2.0*pi*x[2])**2 - sin(2.0*pi*x[2])**2)

class u6Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return sin(2.0*pi*X[0])**2 + cos(2.0*pi*X[1])**2 + sin(2.0*pi*X[2])**2 + X[0]+X[1]+X[2]+5.0
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = [1.0+4.0*pi*cos(2.0*pi*X[0])*sin(2.0*pi*X[0]),
              1.0-4.0*pi*cos(2.0*pi*X[1])*sin(2.0*pi*X[1]),
              1.0+4.0*pi*cos(2.0*pi*X[2])*sin(2.0*pi*X[2])]
        return numpy.reshape(numpy.array(du),(nd,))
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC6(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [1.0]:
        return lambda x,t: u6Ex().uOfXT(x,t)

# def getAdvFluxBC6(x,flag):
#     pass
# def getDiffFluxBC6(x,flag):
#     if x[2] == 0.0:
#         n = numpy.zeros((nd,),'d'); n[2]=-1.0
#         return lambda x,t: numpy.dot(velEx(u6Ex(),a6).uOfXT(x,t),n)
def getDiffFluxBC6(x,flag):
    pass
def getAdvFluxBC6(x,flag):
    if x[2] == 0.0:
        n = numpy.zeros((nd,),'d'); n[2]=-1.0
        return lambda x,t: numpy.dot(velEx(u6Ex(),a6).uOfXT(x,t),n)

##################################################
##################################################
#u = x*x + y*y + z*z, a_00 = x*x + 5, a_11 = y*y + 5.0, a_22 = z*z + 5.0
#f = -(6*x*x + 6*y*y + 6*z*z + 30)
#dirichlet bc's everywhere

def a7(x):
    return numpy.array([[x[0]**2 + 5.0,0.0,0.0],[0.0,x[1]**2 + 5.0,0.0],
                          [0.0,0.0,x[2]**2+5.0]],'d')
def f7(x):
    return -6.0*x[0]**2 -6.0*x[1]**2 - 6.0*x[2]**2 -30.0

class u7Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2+x[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:nd],(nd,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC7(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0]:
        return lambda x,t: u7Ex().uOfXT(x,t)
def getDiffFluxBC7(x,flag):
    pass
def getAdvFluxBC7(x,flag):
    if x[2] == 1.0:
        n = numpy.zeros((nd,),'d'); n[2]=1.0
        return lambda x,t: numpy.dot(velEx(u7Ex(),a7).uOfXT(x,t),n)

##################################################
#u = x*x + y*y + z*z, a_00 = x*x*y + 5, a_11 = y*y*z + 5.0, a_22 = z*z*x + 5.0
#f = -(6*x*x*y + 6*y*y*z + 6*z*z*x + 30)
#dirichlet bc's everywhere

def a8(x):
    return numpy.array([[x[1]*x[0]**2 + 5.0,0.0,0.0],[0.0,x[2]*x[1]**2 + 5.0,0.0],
                          [0.0,0.0,x[0]*x[2]**2+5.0]],'d')
def f8(x):
    return -6.0*x[1]*x[0]**2 -6.0*x[2]*x[1]**2 - 6.0*x[0]*x[2]**2 -30.0

class u8Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2+x[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*numpy.reshape(X[0:nd],(nd,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC8(x,flag):
    if abs(x[0]-0.0) < 1.0e-4 or abs(x[0]-1.0) < 1.0e-4 or abs(x[1]-0.0) < 1.0e-4 or abs(x[1]-1.0) < 1.0e-4 or abs(x[2]-0.0)< 1.0e-4:
        return lambda x,t: u8Ex().uOfXT(x,t)
def getDiffFluxBC8(x,flag):
    pass
def getAdvFluxBC8(x,flag):
    if abs(x[2]-1.0) < 1.0e-4:
        n = numpy.zeros((nd,),'d'); n[2]=1.0
        return lambda x,t: numpy.dot(velEx(u8Ex(),a8).uOfXT(x,t),n)




nc = 1

#ex 0
analyticalSolution = {0:u0Ex()}
dirichletConditions = {0:getDBC0}
aOfX = {0:a0}; fOfX = {0:f0}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC0}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC0}}
#ex 1
#analyticalSolution = {0:u1Ex()}
#dirichletConditions = {0:getDBC1}
#aOfX = {0:a1}; fOfX = {0:f1}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC1}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC1}}
#ex 2
analyticalSolution = {0:u2Ex()}
dirichletConditions = {0:getDBC2}
aOfX = {0:a2}; fOfX = {0:f2}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC2}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC2}}
#ex 3
#analyticalSolution = {0:u3Ex()}
#dirichletConditions = {0:getDBC3}
#aOfX = {0:a3}; fOfX = {0:f3}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC3}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC3}}
#ex 4
#analyticalSolution = {0:u4Ex()}
#dirichletConditions = {0:getDBC4}
#aOfX = {0:a4}; fOfX = {0:f4}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC4}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC4}}
#ex 5
analyticalSolution = {0:u5Ex()}
dirichletConditions = {0:getDBC5}
aOfX = {0:a5}; fOfX = {0:f5}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC5}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC5}}
#ex 6
analyticalSolution = {0:u6Ex()}
dirichletConditions = {0:getDBC6}
aOfX = {0:a6}; fOfX = {0:f6}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC6}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC6}}
#ex 7
analyticalSolution = {0:u7Ex()}
dirichletConditions = {0:getDBC7}
aOfX = {0:a7}; fOfX = {0:f7}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC7}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC7}}
#ex 8
analyticalSolution = {0:u8Ex()}
dirichletConditions = {0:getDBC8}
aOfX = {0:a8}; fOfX = {0:f8}
advectiveFluxBoundaryConditions =  {0:getAdvFluxBC8}
diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC8}}

analyticalSolutionVelocity = {0:velEx(analyticalSolution[0],aOfX[0])}
#mwf hack to test projections
#analyticalSolution = None
#analyticalSolutionVelocity = None

coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,nc,nd)

coefficients.variableNames=['u0']

   

fluxBoundaryConditions = {0:'setFlow'}



