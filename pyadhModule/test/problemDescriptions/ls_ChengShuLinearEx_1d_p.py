from pyadh import *
from pyadh.default_p import *
from math import *

nd =1

def generateAplot(soln,xLeft,xRight,nnx,t,filename='ChengShu_ex2_exactSoln.dat'):
    """
    save exact solution to a file 
    """
    dx = (xRight-xLeft)/(nnx-1.)
    fout = open(filename,'w')
    for i in range(nnx):
        x = (xLeft+dx*i,0.0,0.0)
        u = soln.uOfXT(x,t)
        fout.write('%12.5e  %12.5e \n' % (x[0],u))
    #
    fout.close()

"""
Example 1
phi_t + sin(x)phi_x =0
phi(x,0) = sin(x)
periodic bcs on [0,2pi]
"""
L = [2.0*math.pi,1.0,1.0]
class Example1IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return math.sin(x[0])
class Example1Soln:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return math.sin(2.0*math.atan(math.exp(-t)*math.tan(0.5*x[0])))
T1 = 1.0
class Example1HJ(TransportCoefficients.TC_base):
    def __init__(self):
        mass={0:{0:'linear'}}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={0:{0:'linear'}}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
    def evaluate(self,t,c):
        c[('m',0)][:]=c[('u',0)]
        c[('dm',0,0)][:]=1.0
        c[('H',0)].flat[:]=c[('grad(u)',0)].flat[:]*numpy.sin(c['x'].take([0],axis=-1).flat[:])
        c[('dH',0,0)].flat[:]=numpy.sin(c['x'].take([0],axis=-1).flat)
        #c[('velocity',0)][:]=1.0
        if c.has_key(('velocity',0)):
            c[('velocity',0)][:]=c[('dH',0,0)].flat[:]

"""
Example 2
phi_t + sign(cos(x))phi_x = 0
phi(x,0) = sin(x)
periodic bcs on [0,2pi]

"""
L = [2.0*math.pi,1.0,1.0]
   
class Example2IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return math.sin(x[0])
class Example2Soln:
    def __init__(self):
        #pass
        #mwf having problems getting code to plot exact solution
        generateAplot(self,0.0,2.0*math.pi,101,1.0,filename='ChengShu_ex2_exactSoln.dat')
    def uOfXT(self,x,t):
        if 0. <= t and t <= 0.5*math.pi:
            return self.ut0_pi2(x,t)
        if 0.5*math.pi < t and t <= math.pi:
            return self.utpi2_pi(x,t)
        return -1.0
    def ut0_pi2(self,x,t):
        if 0.0 <= x[0] and x[0] <= 0.5*math.pi:
            return math.sin(x[0]-t)
        if 0.5*math.pi < x[0] and x[0] <= 1.5*math.pi-t:
            return math.sin(x[0]+t)
        if 1.5*math.pi-t < x[0] and x[0] <= 1.5*math.pi+t:
            return -1.0
        if 1.5*math.pi+t < x[0] and x[0] <= 2.0*math.pi+1.0e-6:
            return sin(x[0]-t)
    def utpi2_pi(self,x,t):
        if 0.0 <= x[0] and x[0] <= t - 0.5*math.pi:
            return -1.0
        if t - 0.5*math.pi < x[0] and x[0] <= 0.5*math.pi:
            return math.sin(x[0]-t)
        if 0.5*math.pi < x[0] and x[0] <= 1.5*math.pi-t:
            return math.sin(x[0]+t)
        if 1.5*math.pi-t < x[0] and x[0] <= 2.0*math.pi+1.0e-6:
            return -1.0
         
T2 = 1.0

class Example2HJ(TransportCoefficients.TC_base):
    def __init__(self):
        mass={0:{0:'linear'}}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={0:{0:'linear'}}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
    def evaluate(self,t,c):
        c[('m',0)][:]=c[('u',0)]
        c[('dm',0,0)][:]=1.0
        c[('H',0)].flat[:]=c[('grad(u)',0)].flat[:]*numpy.sign(numpy.cos(c['x'].take([0],axis=-1).flat[:]))
        c[('dH',0,0)].flat[:]=numpy.sign(numpy.cos(c['x'].take([0],axis=-1).flat[:]))
        #c[('velocity',0)][:]=1.0
        if c.has_key(('velocity',0)):
            c[('velocity',0)][:]=c[('dH',0,0)].flat[:]


"""
Example 3
phi_t + (phi_x)^2/2 = 0
phi(x,0) = sin(x)
periodic bcs on [0,2pi]

"""
L = [2.0*math.pi,1.0,1.0]
class Example3IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return math.sin(x[0])
#need to compute actual solution
Example3Soln = Example3IC

#still smooth        
#T3 = 0.5
T3 = 1.5

class Example3HJ(TransportCoefficients.TC_base):
    def __init__(self):
        mass={0:{0:'linear'}}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={0:{0:'nonlinear'}}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
    def evaluate(self,t,c):
        c[('m',0)][:]=c[('u',0)]
        c[('dm',0,0)][:]=1.0
        c[('H',0)].flat[:]=0.5*c[('grad(u)',0)].flat[:]*c[('grad(u)',0)].flat[:]
        c[('dH',0,0)].flat[:]=c[('grad(u)',0)].flat[:]
        #c[('velocity',0)][:]=1.0
        if c.has_key(('velocity',0)):
            c[('velocity',0)][:]=c[('dH',0,0)].flat[:]


"""
Example 4
phi_t + (phi_x)^2/2 = 0
phi(x,0) = pi - x   0 <= x <= pi
           x  - pi  pi <= x <= 2pi
periodic bcs on [0,2pi]

"""
L = [2.0*math.pi,1.0,1.0]
class Example4IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if 0.0 <= x[0] and x[0] <= math.pi:
            return math.pi - x[0]
        else:
            return x[0] - math.pi
        
#need to compute actual solution
Example4Soln = Example4IC

#still smooth        
T4 = 1.

Example4HJ = Example3HJ

"""
Example 5
phi_t + |phi_x| = 0
phi(x,0) = sin(x)
periodic bcs on [0,2pi]

"""
L = [2.0*math.pi,1.0,1.0]
class Example3IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return math.sin(x[0])
#need to compute actual solution
Example5Soln = Example2Soln

class Example5HJ(TransportCoefficients.TC_base):
    def __init__(self):
        mass={0:{0:'linear'}}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={0:{0:'nonlinear'}}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
    def evaluate(self,t,c):
        c[('m',0)][:]=c[('u',0)]
        c[('dm',0,0)][:]=1.0
        c[('H',0)].flat[:]=numpy.sqrt(c[('grad(u)',0)].flat[:]*c[('grad(u)',0)].flat[:])
        c[('dH',0,0)].flat[:]=c[('grad(u)',0)].flat[:]/(numpy.abs(c[('grad(u)',0)].flat[:])+1.0e-12)
        if c.has_key(('velocity',0)):
            c[('velocity',0)][:]=c[('dH',0,0)].flat[:]

T5 = 1.0


#analyticalSolution = {0:Example1Soln()}
#coefficients = Example1HJ()
#T = T1

analyticalSolution = {0:Example2Soln()}
coefficients = Example2HJ()
T = T2

#analyticalSolution = {0:Example3Soln()}
#coefficients = Example3HJ()
#T = T3

#analyticalSolution = {0:Example4Soln()}
#coefficients = Example4HJ()
#T = T4

#analyticalSolution = {0:Example5Soln()}
#coefficients = Example5HJ()
#T = T5

name = 'ls_ChengShu_ex4_dgp2_Jaffre_nuc0_1_nn41'

coefficients.variableNames=['phi']

#now define the Dirichlet boundary conditions

def getDBC(x,tag):
    pass
    
dirichletConditions = {0:getDBC}

def getPDBC(x,tag):
    if x[0] == 0.0 or x[0] == L[0]:
        return numpy.array([0.0,0.0,0.0])

#periodic boundary conditions

periodicDirichletConditions = {0:getPDBC}

initialConditions  = {0:analyticalSolution[0]}

def getAFBC(x,tag):
    pass

advectiveFluxBoundaryConditions =  {}#can't have for dg {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}

## @}
