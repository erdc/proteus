from pyadh import *
from pyadh.default_p import *
"""
Heterogeneous Poisson's equations for one component (uncoupled) in 2D
"""

##\page Tests Test Problems 
# \ref poisson_1c_het_2d_p.py "Heterogeneous Poisson's equation for 1 components"
#

##\ingroup test
#\file poisson_1c_het_2d_p.py
#
#\brief Heterogensou Poisson's equations for one component (uncoupled) in 2D

nd = 2

initialConditions = None
sd=True#False
Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0

class velEx:
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = numpy.reshape(self.aex(X),(2,2))
        #mwf debug
        #print """du =%s a= %s v= %s """ % (du,A,-numpy.dot(A,du))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)


##################################################
#u = 1-x_{i}, a = I, f = 0 
#dir bc's along x, no flow along other axis
ax = 1 #0 1
ay =  (ax+1) % 2
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
        du = numpy.zeros((2,),'d')
        du[ax] = -1.0
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)
    
def getDBC0(x,flag):
    if x[ax] == 0.0:
        return lambda x,t: 1.0
    if x[ax] == 1.0:
        return lambda x,t: 0.0
#     if x[ax+1] == 0.0:
#         return lambda x,t: 1.0-x[0]
#     if x[ax+1] == 1.0:
#         return lambda x,t: 1.0-x[0]

def getAdvFluxBC0(x,flag):
    if x[ay] == 0.0:
        return lambda x,t: 0.0
    if x[ay] == 1.0:
        return lambda x,t: 0.0

def getDiffFluxBC0(x,flag):
    if x[ay] == 0.0:
        return lambda x,t: 0.0
    if x[ay] == 1.0:
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
        du = numpy.ones((2,),'d')
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC1(x,flag):
    if x[0] in [1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u1Ex().uOfXT(x,t)

def getDiffFluxBC1(x,flag):
   if x[0] == 0.0:
       n = numpy.array([-1.0,0.0],'d')
       #mwf debug
       #print """getAdvFlux1 x= %s q.n= %s """ % (x,numpy.dot(velEx(u1Ex(),a1).uOfXT(x,0.0),n))
       return lambda x,t: numpy.dot(velEx(u1Ex(),a1).uOfXT(x,t),n)

def getAdvFluxBC1(x,flag):
    pass

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
#u = sin^2(2\pi x) + cos^2(2\pi y) + x+y + 5.0, a = I,
#f = -8\pi^2cos^2(2\pi x) + 8\pi^2cos^2(2\pi y)
#    +8\pi^2sin^2(2\pi x) - 8\pi^2sin^2(2\pi y)
#dirichlet bc's everywhere
#cek hack, I've messed this test problem up
def a3(x):
    return Ident.flat[:]
def f3(x):
    #return 8.0*pi**2*(-cos(2.*pi*x[0])**2 + cos(2.0*pi*x[1])**2 + sin(2.0*pi*x[0])**2 - sin(2.0*pi*x[1])**2)
    return 4.0*pi**2 * (sin(2.0*pi*x[0]) + sin(2.0*pi*x[1]))

class u3Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        #return sin(2.0*pi*X[0])**2 + cos(2.0*pi*X[1])**2 + X[0]+X[1]+5.0
        return sin(2.0*pi*X[0]) + sin(2.0*pi*X[1])
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = [1.0+4.0*pi*cos(2.0*pi*X[0])*sin(2.0*pi*X[0]),
              1.0-4.0*pi*cos(2.0*pi*X[1])*sin(2.0*pi*X[1])]
        return numpy.reshape(numpy.array(du),(2,))
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC3(x,flag):
    #if x[0] in [0.0,1.0] or x[1] in [0.0]:
    #    return lambda x,t: u3Ex().uOfXT(x,t)
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u3Ex().uOfXT(x,t)
def getDiffFluxBC3(x,flag):
    pass
def getAdvFluxBC3(x,flag):
    pass
#     if x[1] == 1.0:
#         n = numpy.zeros((nd,),'d'); n[1]=1.0
#         return lambda x,t: numpy.dot(velEx(u3Ex(),a3).uOfXT(x,t),n)

##################################################
##################################################
#u = x + y, a_00 = sin^2(2\pi x) + 5, a_11 = cos^2(2\pi y) + 5.0
#f = 4\pi*(-cos(2\pi*x)sin(2\pi*x) + cos(2\pi*y)sin(2\pi*y))
#dirichlet bc's everywhere

def a4(x):
    return numpy.array([[sin(2.0*pi*x[0])**2 + 5.0,0.0],[0.0,cos(2.0*pi*x[1])**2 + 5.0]],'d')

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
        du = numpy.ones((2,),'d')
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC4(x,flag):
    if x[0] in [0.0,1.0] or x[1] in [1.0]:
        return lambda x,t: u4Ex().uOfXT(x,t)
#for pwl might be easier to enforce boundary flux as advective flux
def getAdvFluxBC4(x,flag):
    pass
def getDiffFluxBC4(x,flag):
    if x[1] == 0.0:
        n = numpy.zeros((nd,),'d'); n[1]=-1.0
        return lambda x,t: numpy.dot(velEx(u4Ex(),a4).uOfXT(x,t),n)
#def getAdvFluxBC4(x,flag):
#    pass
#def getDiffFluxBC4(x,flag):
#     if x[1] == 0.0:
#         n = numpy.zeros((nd,),'d'); n[1]=-1.0
#         return lambda x,t: numpy.dot(velEx(u4Ex(),a4).uOfXT(x,t),n)
    
##################################################
##################################################
#u = x*x + y*y, a_00 = x + 5, a_11 = y + 5.0
#f = -2*x -2*(5+x) -2*y-2*(5+y)
#

def a5(x):
    return numpy.array([[x[0] + 5.0,0.0],[0.0,x[1] + 5.0]],'d')
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
        du = 2.0*numpy.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC5(x,flag):
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u5Ex().uOfXT(x,t)
def getAdvFluxBC5(x,flag):
   pass
def getDiffFluxBC5(x,flag):
    if x[0] == 1.0:
        n = numpy.zeros((nd,),'d'); n[0]=1.0
        return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)


#def getAdvFluxBC5(x,flag):
#    pass
#def getDiffFluxBC5(x,flag):
#    if x[0] == 1.0:
#        n = numpy.zeros((nd,),'d'); n[0]=1.0
#        return lambda x,t: numpy.dot(velEx(u5Ex(),a5).uOfXT(x,t),n)
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
    return numpy.array([[x[0]**2 + 5.0,0.0],[0.0,x[1]**2 + 5.0]],'d')
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
        return numpy.reshape(numpy.array(du),(2,))
    def duOfXT(self,X,T):
        return self.duOfX(X)
eps=1.0e-8
def getDBC6(x,flag):
    #if x[0] in [0.0,1.0] or x[1] in [0.0]:
    if abs(x[0]-0.0) < eps or abs(x[0]-1.0) < eps or abs(x[1]-0.0) < eps:
        return lambda x,t: u6Ex().uOfXT(x,t)
def getAdvFluxBC6(x,flag):
    if abs(x[1]-1.0) < eps:
        n = numpy.zeros((nd,),'d'); n[1]=1.0
        return lambda x,t: numpy.dot(velEx(u6Ex(),a6).uOfXT(x,t),n)
def getDiffFluxBC6(x,flag):
    if abs(x[1] - 1.0) < eps:
        return lambda x,t: 0.0

#def getDiffFluxBC6(x,flag):
#    pass
# def getAdvFluxBC6(x,flag):
#     pass
# def getDiffFluxBC6(x,flag):
#     if x[1] == 1.0:
#         n = numpy.zeros((nd,),'d'); n[1]=1.0
#         return lambda x,t: numpy.dot(velEx(u6Ex(),a6).uOfXT(x,t),n)
##################################################
##################################################
#u = x*x + y*y, a_00 = x*x + 5, a_11 = y*y + 5.0
#f = -(6*x*x + 6*y*y + 20)
#

def a7(x):
    return numpy.array([[x[0]**2 + 5.0,0.0],[0.0,x[1]**2 + 5.0]],'d')
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
        du = 2.0*numpy.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC7(x,flag):
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u7Ex().uOfXT(x,t)
def getAdvFluxBC7(x,flag):
    if x[0] == 1.0:
        n = numpy.zeros((nd,),'d'); n[0]=1.0
        return lambda x,t: numpy.dot(velEx(u7Ex(),a7).uOfXT(x,t),n)
def getDiffFluxBC7(x,flag):
    pass
##################################################
##################################################
#u = x*x + y*y, a_00 = x*x*y + 5, a_11 = y*y*x + 5.0
#f = -(6*x*x*y + 6*y*y*x + 20)
#

def a8(x):
    return numpy.array([[x[1]*x[0]**2 + 5.0,0.0],[0.0,x[0]*x[1]**2 + 5.0]],'d')
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
        du = 2.0*numpy.reshape(X[0:2],(2,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC8(x,flag):
    if x[0] in [0.0] or x[1] in [0.0,1.0]:
        return lambda x,t: u8Ex().uOfXT(x,t)
def getAdvFluxBC8(x,flag):
    if x[0] == 1.0:
        n = numpy.zeros((nd,),'d'); n[0]=1.0
        return lambda x,t: numpy.dot(velEx(u8Ex(),a8).uOfXT(x,t),n)
def getDiffFluxBC8(x,flag):
    pass

# def getDBC8(x,flag):
#     if x[0] in [0.0,1.0] or x[1] in [0.0,1.0]:
#         return lambda x,t: u8Ex().uOfXT(x,t)
# def getAdvFluxBC8(x,flag):
#     pass
# def getDiffFluxBC8(x,flag):
#     pass
##################################################

nc = 1

#ex 0
#analyticalSolution = {0:u0Ex()}
#dirichletConditions = {0:getDBC0}
#aOfX = {0:a0}; fOfX = {0:f0}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC0}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC0}}
#
#
#ex 1
# analyticalSolution = {0:u1Ex()}
# dirichletConditions = {0:getDBC1}
# aOfX = {0:a1}; fOfX = {0:f1}
# advectiveFluxBoundaryConditions =  {0:getAdvFluxBC1}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC1}}
#
#ex 2
#analyticalSolution = {0:u2Ex()}
#dirichletConditions = {0:getDBC2}
#aOfX = {0:a2}; fOfX = {0:f2}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC2}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC2}}
#
#ex 3
# analyticalSolution = {0:u3Ex()}
# dirichletConditions = {0:getDBC3}
# aOfX = {0:a3}; fOfX = {0:f3}
# advectiveFluxBoundaryConditions =  {0:getAdvFluxBC3}
# diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC3}}

#ex 4
#analyticalSolution = {0:u4Ex()}
#dirichletConditions = {0:getDBC4}
#aOfX = {0:a4}; fOfX = {0:f4}
#advectiveFluxBoundaryConditions =  {0:getAdvFluxBC4}
#diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC4}}
#
#ex 5
# analyticalSolution = {0:u5Ex()}
# dirichletConditions = {0:getDBC5}
# aOfX = {0:a5}; fOfX = {0:f5}
# advectiveFluxBoundaryConditions =  {0:getAdvFluxBC5}
# diffusiveFluxBoundaryConditions = {0:{0:getDiffFluxBC5}}
#
#ex 6
analyticalSolution = dict((i,u6Ex()) for i in range(nc))   #{0:u6Ex(),1:u6Ex()}
dirichletConditions = dict((i,getDBC6) for i in range(nc)) #{0:getDBC6,1:getDBC6}
aOfX = dict((i,a6) for i in range(nc)) ; fOfX = dict((i,f6) for i in range(nc)) #{0:a6,1:a6}; fOfX = {0:f6,1:f6}
advectiveFluxBoundaryConditions =  dict((i,getAdvFluxBC6) for i in range(nc))   #{0:getAdvFluxBC6,1:getAdvFluxBC6}
diffusiveFluxBoundaryConditions = dict((i,{i:getDiffFluxBC6}) for i in range(nc))#{0:{0:getDiffFluxBC6},1:{1:getDiffFluxBC6}}

#ex 7
#analyticalSolution = dict((i,u7Ex()) for i in range(nc))   #{0:u7Ex(),1:u7Ex()}
#dirichletConditions = dict((i,getDBC7) for i in range(nc)) #{0:getDBC7,1:getDBC7}
#aOfX = dict((i,a7) for i in range(nc)) ; fOfX = dict((i,f7) for i in range(nc)) #{0:a7,1:a7}; fOfX = {0:f7,1:f7}
#advectiveFluxBoundaryConditions =  dict((i,getAdvFluxBC7) for i in range(nc))   #{0:getAdvFluxBC7,1:getAdvFluxBC7}
#diffusiveFluxBoundaryConditions = dict((i,{i:getDiffFluxBC7}) for i in range(nc))#{0:{0:getDiffFluxBC7},1:{1:getDiffFluxBC7}}


analyticalSolutionVelocity = dict((i,velEx(analyticalSolution[i],aOfX[i])) for i in range(nc))#{0:velEx(analyticalSolution[0],aOfX[0]),1:velEx(analyticalSolution[1],aOfX[1])}
#mwf hack to test projections
#analyticalSolutionVelocity = None
#analyticalSolution = None
coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,nc,nd)

coefficients.variableNames=list('u%d' % i for i in range(nc))#['u0','u1']

   
fluxBoundaryConditions = dict((i,'setFlow') for i in range(nc))#{0:'setFlow',1:'setFlow'}



