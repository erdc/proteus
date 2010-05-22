from pyadh import *
from pyadh.default_p import *
"""
linear diffusion-reaction system with two-components 2D.
equations are uncoupled. includes 'constant' gravity type term 
"""

## \page Tests Test Problems
#\ref dr_2c_di_het_ss_3d_p.py "Linear diffusion-reaction, two-component, diagonal system (3D)"

##\ingroup test
#\file dr_2c_di_het_ss_3d_p.py
#
#\brief linear diffusion-reaction system with two-components in 3D. includes
#  constant 'gravity-like' term to test models for linear flow 
#\todo finish dr_2c_di_het_ss_3d_p.py impl,doc
nd = 3

initialConditions = None

#identity tensor in 3d
Ident = Numeric.zeros((nd,nd),Numeric.Float)
Ident[0,0]=1.0; Ident[1,1] = 1.0; Ident[2,2] = 1.0;
#direction of gravity in system
gravityU = Numeric.zeros((nd,),Numeric.Float)
gravityU[-1] = -1.0

class velEx:
    def __init__(self,duex,aex,bex):
        self.duex = duex
        self.aex = aex
        self.bex = bex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = self.aex(X)
        b  = self.bex(X)
        #mwf debug
        #print """du =%s a= %s b=%s  v= %s """ % (du,A,b,-Numeric.dot(A,(du-b)))
        return -Numeric.dot(A,(du-b))
    def uOfXT(self,X,T):
        return self.uOfX(X)

##################################################
#u = 1- x_{i}, a = I, b = [0,0,-1], c = 1.0 f = 1- x_i 
#dir bc's everywhere
ax = 0
def a0(x):
    return Ident
def b0(x):
    return gravityU
def c0(x):
    return 1.0
def f0(x):
    return 1.-x[ax]


class u0Ex:
    def __init__(self):
        self.ax = ax
    def uOfX(self,X):
        return 1.-X[self.ax]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.zeros((nd,),Numeric.Float)
        du[self.ax]=-1.0
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)
    
def getDBC0(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u0Ex().uOfXT(x,t)




##################################################
#u = \sum x_{i}, a = I, b = [0,0,-1], c = 1.0 f = \sum x_i 
#dir bc's everywhere
def a1(x):
    return Ident
def b1(x):
    return gravityU
def c1(x):
    return 1.0
def f1(x):
    return x[0]+x[1]+x[2]


class u1Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return X[0]+X[1]+X[2]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.ones((nd,),Numeric.Float)
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)
    
def getDBC1(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u1Ex().uOfXT(x,t)

##################################################
#u = \sum x_{i}^2, a = I, f = -2*d 
#dirichlet bc's everywhere
def a2(x):
    return Ident
def b2(x):
    return gravityU
def c2(x):
    return 1.0
def f2(x):
    return x[0]**2+x[1]**2+x[2]**2 -2.0*nd


class u2Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return X[0]**2+X[1]**2+X[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*X[0:nd]
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)
    
def getDBC2(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u2Ex().uOfXT(x,t)

##################################################
#u = sin^2(2\pi x) + cos^2(2\pi y) + sin^2(2\pi z) + x+y+z + 5.0, a = I, b = [0,0,-1], c = 1.0
#f = -8\pi^2cos^2(2\pi x) + 8\pi^2cos^2(2\pi y) -8\pi^2cos^2(2\pi z)
#    +8\pi^2sin^2(2\pi x) - 8\pi^2sin^2(2\pi y) +8\pi^2sin^2(2\pi z)
#    + sin^2(2\pi x) + cos^2(2\pi y) + sin^2(2\pi z) + x+y+z + 5.0
def a3(x):
    return Ident
def b3(x):
    return gravityU
def c3(x):
    return 1.0
def f3(x):
    return 8.0*pi**2*(-cos(2.*pi*x[0])**2 + cos(2.0*pi*x[1])**2 + sin(2.0*pi*x[0])**2 - sin(2.0*pi*x[1])**2 -cos(2.*pi*x[2])**2 + sin(2.0*pi*x[2])**2) + sin(2.0*pi*x[0])**2 + cos(2.0*pi*x[1])**2 + sin(2.0*pi*x[2]**2) + x[0]+x[1]+x[2]+5.0


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
        return Numeric.reshape(Numeric.array(du),(nd,))
    def duOfXT(self,X,T):
        return self.duOfX(X)
    
def getDBC3(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u3Ex().uOfXT(x,t)


########################################################################
#u = x + y + z, a_00 = sin^2(2\pi x) + 5, a_11 = cos^2(2\pi y) + 5.0 
# a_22 = sin^2(2\pi z) + 5  c = y+1
#f = 4\pi*(-cos(2\pi*x)sin(2\pi*x) + cos(2\pi*y)sin(2\pi*y)
#          -cos(2\pi*z)sin(2\pi*z))
#   - 4\pi*sin(2\pi*z)*cos(2\pi*z) + (x+y+z)(y+1) 
#dirichlet bc's everywhere

def a4(x):
    return Numeric.array([[sin(2.0*pi*x[0])**2 + 5.0,0.0,0.0],
                          [0.0,cos(2.0*pi*x[1])**2 + 5.0,0.0],
                          [0.0,0.0,sin(2.0*pi*x[2])**2 + 5]],Numeric.Float)
def b4(x):
    return gravityU
def c4(x):
    return 1.0+x[1]
def f4(x):
    return 4.0*pi*(-cos(2.0*pi*x[0])*sin(2.0*pi*x[0]) + \
                    cos(2.0*pi*x[1])*sin(2.0*pi*x[1])  - \
                    cos(2.0*pi*x[2])*sin(2.0*pi*x[2])) - \
                   4.0*pi*cos(2.0*pi*x[2])*sin(2.0*pi*x[2]) + (x[0]+x[1]+x[2])*(x[1]+1.0)

class u4Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]+x[1]+x[2]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.ones((nd,),Numeric.Float)
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC4(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u4Ex().uOfXT(x,t)
#
########################################################################
#u = x + y + z, a_00 = sin^2(2\pi x) + 5, a_11 = cos^2(2\pi y) + 5.0
#    a_22 = sin^2(2\pi z) + 5
#    b = [0,0,-z], c =1
#f = 4\pi*(-cos(2\pi*x)sin(2\pi*x) + cos(2\pi*y)sin(2\pi*y)
#          -cos(2\pi*z)sin(2\pi*z))
#   - 4\pi*z*sin(2\pi*z)*cos(2\pi*z) - sin^2(2\pi z) -5 + x+y+z
#dirichlet bc's everywhere

def a4a(x):
    return Numeric.array([[sin(2.0*pi*x[0])**2 + 5.0,0.0,0.0],
                          [0.0,cos(2.0*pi*x[1])**2 + 5.0,0.0],
                          [0.0,0.0,sin(2.0*pi*x[2])**2 + 5]],Numeric.Float)
def b4a(x):
    b = Numeric.zeros((nd,),Numeric.Float)
    b[2]=-x[2]
    return b
def c4a(x):
    return 1.0
def f4a(x):
    return 4.0*pi*(-cos(2.0*pi*x[0])*sin(2.0*pi*x[0]) + \
                    cos(2.0*pi*x[1])*sin(2.0*pi*x[1])  - \
                    cos(2.0*pi*x[2])*sin(2.0*pi*x[2])) - \
                   4.0*pi*x[2]*cos(2.0*pi*x[2])*sin(2.0*pi*x[2]) - \
                   sin(2.0*pi*x[2])**2 - 5.0 + x[0] + x[1] + x[2]


class u4aEx:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]+x[1]+x[2]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.ones((nd,),Numeric.Float)
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC4a(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u4aEx().uOfXT(x,t)
#
########################################################################
#u = x + y, a_00 = sin^2(2\pi x) + 5, a_11 = cos^2(2\pi y) + 5.0,
#    a_22 = sin^2(2\pi z) + 5
#    b = [0,0,-z], c = y
#f = 4\pi*(-cos(2\pi*x)sin(2\pi*x) + cos(2\pi*y)sin(2\pi*y)
#          -cos(2\pi*z)sin(2\pi*z))
#   - 4\pi*z*sin(2\pi*z)*cos(2\pi*z) - sin^2(2\pi z) -5 + (x+y+z)*x
#dirichlet bc's everywhere

def a4b(x):
    return Numeric.array([[sin(2.0*pi*x[0])**2 + 5.0,0.0,0.0],
                          [0.0,cos(2.0*pi*x[1])**2 + 5.0,0.0],
                          [0.0,0.0,sin(2.0*pi*x[2])**2 + 5]],Numeric.Float)

def b4b(x):
    b = Numeric.zeros((nd,),Numeric.Float)
    b[2]=-x[2]
    return b
def c4b(x):
    return x[1]
def f4b(x):
    return 4.0*pi*(-cos(2.0*pi*x[0])*sin(2.0*pi*x[0]) + \
                    cos(2.0*pi*x[1])*sin(2.0*pi*x[1])  - \
                    cos(2.0*pi*x[2])*sin(2.0*pi*x[2])) - \
                   4.0*pi*x[2]*cos(2.0*pi*x[2])*sin(2.0*pi*x[2]) - \
                   sin(2.0*pi*x[2])**2 - 5.0 + (x[0]+x[1]+x[2])*x[0]


class u4bEx:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]+x[1]+x[2]
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = Numeric.ones((nd,),Numeric.Float)
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC4b(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u4bEx().uOfXT(x,t)
#
##################################################
#u = x*x + y*y + z*z, a_00 = x + 5, a_11 = y + 5.0, a_22 = z + 5.0
#    b = [0,0,-z], c = (y+1)
#f = -2*x -2*(5+x) -2*y-2*(5+y) -2*z -2*(5+z) -2*z - 5 + (y+1)*(x*x+y*y+z*z)
#dirichlet bc's everywhere
def a5(x):
    return Numeric.array([[x[0] + 5.0,0.0,0.0],
                          [0.0,x[1] + 5.0,0.0],
                          [0.0,0.0,x[2] + 5.0]],Numeric.Float)
def b5(x):
    b=Numeric.zeros((nd,),Numeric.Float)
    b[2] = -x[2]
    return b
def c5(x):
    return x[1]+1.0
def f5(x):
    return -2.0*x[0] -2*(5.+x[0]) -2.*x[1]-2.*(5.+x[1]) -2.0*x[2] -2*(5.+x[2]) \
           -2.0*x[2]-5.0 +(x[1]+1.)*(x[0]**2+x[1]**2+x[2]**2)

class u5Ex:
    def __init__(self):
        pass
    def uOfX(self,x):
        return x[0]**2+x[1]**2+x[2]**2
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = 2.0*Numeric.reshape(X[0:nd],(nd,))
        return du
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC5(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u5Ex().uOfXT(x,t)


##################################################
#u = sin^2(2\pi x) + cos^2(2\pi y) + sin^2(2\pi z) + x+y+z+ 5.0,
#a_00 = x*x + 5.0 , a_11 = y*y + 5.0, a_22 = z*z + 5.0
#b = [0,0,-z], c = x 
#f = -2.*x*(1. + 4*pi*Cos(2*pi*x)*Sin(2*pi*x)) - 
#   (5 + x*x)*(8*pi*pi*Cos^2(2*pi*x) - 8*pi*pi*Sin^2(2*pi*x)) - 
#   2.*y*(1. - 4*pi*Cos(2*pi*y)*Sin(2*pi*y)) - 
#   (5 + 1.*y*y)*(-8*pi*pi*Cos^2(2*pi*y) + 8*pi*pi*Sin^2(2*pi*y)) -
#   2.*z*(1. + 4*pi*Cos(2*pi*z)*Sin(2*pi*z)) - 
#   (5 + z*z)*(8*pi*pi*Cos^2(2*pi*z) - 8*pi*pi*Sin^2(2*pi*z))
#   -3*z*z - 5 + x*(sin^2(2*pi*x) + cos^2(2*pi*y) + sin^2(2\pi z) + x + y + z + 5.)
########################################################################
#dirichlet bc's everywhere

def a6(x):
    return Numeric.array([[x[0]**2 + 5.0,0.0,0.0],
                          [0.0,x[1]**2 + 5.0,0.0],
                          [0.0,0.0,x[2]**2 + 5.0]],Numeric.Float)
def b6(x):
    b=Numeric.zeros((nd,),Numeric.Float)
    b[2] = -x[2]
    return b
def c6(x):
    return x[0]
def f6(x):
    return -2.*x[0]*(1. + 4.*pi*cos(2.*pi*x[0])*sin(2.*pi*x[0])) - \
           (5.0 + x[0]**2)*8.*pi*pi*(cos(2.0*pi*x[0])**2 - sin(2.0*pi*x[0])**2) - \
           2.*x[1]*(1. - 4.*pi*cos(2.*pi*x[1])*sin(2.*pi*x[1]))  - \
           (5. + x[1]**2)*8.*pi*pi*(-cos(2.0*pi*x[1])**2 + sin(2.*pi*x[1])**2) - \
           2.*x[2]*(1. + 4.*pi*cos(2.*pi*x[2])*sin(2.*pi*x[2])) - \
           (5.0 + x[2]**2)*8.*pi*pi*(cos(2.0*pi*x[2])**2 - sin(2.0*pi*x[2])**2) - \
           3.0*x[2]**2 - 5.0 + x[0]*(sin(2.0*pi*x[0])**2 + cos(2.*pi*x[1])**2 + \
                                     sin(2.0*pi*x[2])**2 + x[0]+x[1]+x[2]+5.)

class u6Ex:
    def __init__(self):
        pass
    def uOfX(self,X):
        return sin(2.0*pi*X[0])**2 + cos(2.0*pi*X[1])**2 + sin(2.0*pi*X[2])**2 + \
               X[0]+X[1]+X[2]+5.0
    def uOfXT(self,X,T):
        return self.uOfX(X)
    #mwf added for velocity calcs
    def duOfX(self,X):
        du = [1.0+4.0*pi*cos(2.0*pi*X[0])*sin(2.0*pi*X[0]),
              1.0-4.0*pi*cos(2.0*pi*X[1])*sin(2.0*pi*X[1]),
              1.0+4.0*pi*cos(2.0*pi*X[2])*sin(2.0*pi*X[2])]
        return Numeric.reshape(Numeric.array(du),(nd,))
    def duOfXT(self,X,T):
        return self.duOfX(X)

def getDBC6(x):
    if x[0] in [0.0,1.0] or x[1] in [0.0,1.0] or x[2] in [0.0,1.0]:
        return lambda x,t: u6Ex().uOfXT(x,t)




##################################################
class EquationCoefficients(TransportCoefficients.TC_base):
    def __init__(self,aOfX,bOfX,cOfX,fOfX,nc=1):
        self.aOfX = aOfX
        self.bOfX = bOfX
        self.cOfX = cOfX
        self.fOfX = fOfX
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        for i in range(nc):
            diffusion[i] = {i : {i:'constant'}}
            advection[i] = {i : 'constant'} #now include for gravity type terms
            reaction[i]  = {i : 'linear'}
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
            for i in range(len(c[('u',ci)].flat)):
                x = c['x'].flat[3*i:3*(i+1)]
                c[('r',ci)].flat[i] = self.cOfX[ci](x)*c[('u',ci)].flat[i]-self.fOfX[ci](x)
                c[('dr',ci,ci)].flat[i] = self.cOfX[ci](x)
                c[('a',ci,ci)].flat[nd*nd*i:nd*nd*(i+1)] = self.aOfX[ci](x).flat
                c[('f',ci)].flat[nd*i:nd*(i+1)] = Numeric.dot(self.aOfX[ci](x),self.bOfX[ci](x)).flat
            #end i
        #end ci
    #end def

######################################################################
#ex 0
#analyticalSolution = {0:u0Ex(),1:u0Ex()}
#dirichletConditions = {0:getDBC0,1:getDBC0}
#aOfX = {0:a0,1:a0}; bOfX= {0:b0,1:b0}; cOfX = {0:c0,1:c0}; fOfX = {0:f0,1:f0}
#ex 1
#analyticalSolution = {0:u1Ex(),1:u1Ex()}
#dirichletConditions = {0:getDBC1,1:getDBC1}
#aOfX = {0:a1,1:a1}; bOfX= {0:b1,1:b1}; cOfX = {0:c1,1:c1}; fOfX = {0:f1,1:f1}
#ex 2
#analyticalSolution = {0:u2Ex(),1:u2Ex()}
#dirichletConditions = {0:getDBC2,1:getDBC2}
#aOfX = {0:a2,1:a2}; bOfX= {0:b2,1:b2}; cOfX = {0:c2,1:c2}; fOfX = {0:f2,1:f2}
#ex 3
#analyticalSolution = {0:u3Ex(),1:u3Ex()}
#dirichletConditions = {0:getDBC3,1:getDBC3}
#aOfX = {0:a3,1:a3}; bOfX= {0:b3,1:b3}; cOfX = {0:c3,1:c3}; fOfX = {0:f3,1:f3}
#ex 4
#analyticalSolution = {0:u4Ex(),1:u4Ex()}
#dirichletConditions = {0:getDBC4,1:getDBC4}
#aOfX = {0:a4,1:a4}; bOfX= {0:b4,1:b4}; cOfX = {0:c4,1:c4}; fOfX = {0:f4,1:f4}
#ex 4a
#analyticalSolution = {0:u4aEx(),1:u4aEx()}
#dirichletConditions = {0:getDBC4a,1:getDBC4a}
#aOfX = {0:a4a,1:a4a}; bOfX= {0:b4a,1:b4a}; cOfX = {0:c4a,1:c4a}; fOfX = {0:f4a,1:f4a}
#ex 4b
#analyticalSolution = {0:u4bEx(),1:u4bEx()}
#dirichletConditions = {0:getDBC4b,1:getDBC4b}
#aOfX = {0:a4b,1:a4b}; bOfX= {0:b4b,1:b4b}; cOfX = {0:c4b,1:c4b}; fOfX = {0:f4b,1:f4b}
#ex 5
#analyticalSolution = {0:u5Ex(),1:u5Ex()}
#dirichletConditions = {0:getDBC5,1:getDBC5}
#aOfX = {0:a5,1:a5}; bOfX= {0:b5,1:b5}; cOfX = {0:c5,1:c5}; fOfX = {0:f5,1:f5}
#ex 6
analyticalSolution = {0:u6Ex(),1:u6Ex()}
dirichletConditions = {0:getDBC6,1:getDBC6}
aOfX = {0:a6,1:a6}; bOfX= {0:b6,1:b6}; cOfX = {0:c6,1:c6}; fOfX = {0:f6,1:f6}

######################################################################
nc = 2
analyticalSolutionVelocity = {0:velEx(analyticalSolution[0],aOfX[0],bOfX[0]),
                              1:velEx(analyticalSolution[1],aOfX[1],bOfX[1])}



coefficients = EquationCoefficients(aOfX,bOfX,cOfX,fOfX,nc)
coefficients.variableNames=['u0','u1']

   

fluxBoundaryConditions = {0:'noFlow',1:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}


