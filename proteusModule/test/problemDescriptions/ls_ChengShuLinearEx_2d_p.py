from pyadh import *
from pyadh.default_p import *
from math import *

nd = 2
"""
Example 1
phi_t -(y-1) phi_x + (x-1) phi_y = 0
phi(x,0) = 0       0.3 <= r
           0.3-r   0.1 < r < 0.3
           0.2     r <= 0.1
r = sqrt((x-1.4)^2 + (y-1.4)^2)

phi(x,t) =phi_0((x-1) cos(t) + (y-1) sin(t),-(x-1) sin(t) + (y-1) cos(t))

periodic bcs on [0,2]^2
"""

L1 = [2.0,2.0,1.0]
class Example1IC:
    def __init__(self,center):
        self.xC = center[0]; self.yC = center[1]
    def uOfXT(self,x,t):
        r = sqrt((x[0]-self.xC)**2 + (x[1]-self.yC)**2)
        if 0.3 <= r:
            return 0.0
        if 0.1 < r and r < 0.3:
            return 0.3-r
        if r <= 0.1:
            return 0.2
        
class Example1Soln:
    def __init__(self,u0):
        self.u0 = u0
    def uOfXT(self,x,t):
        xt = numpy.zeros((3,),'d')
        xt[0] = (x[0]-1.0)*math.cos(t) + (x[1]-1.)*math.sin(t)
        xt[1] = -(x[0]-1.0)*math.sin(t) + (x[1]-1.)*math.cos(t)
        return self.u0.uOfXT(xt,t)
    
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
        c[('dH',0,0)][:] = numpy.concatenate((-(c['x'].take([1],axis=-1)-1.0),
                                               c['x'].take([0],axis=-1,)-1.0),axis=-1)
        #for i in range(len(c[('u',0)].flat)):
        #    c[('dH',0,0)].flat[nd*i+0] = -(c['x'].flat[3*i+1]-1.0)
        #    c[('dH',0,0)].flat[nd*i+1] =  (c['x'].flat[3*i+0]-1.0)
            
                     
        c[('H',0)].flat[:]=c[('grad(u)',0)].take([1],axis=-1).flat[:]*(c['x'].take([0],axis=-1).flat[:]-1.0) - \
                            c[('grad(u)',0)].take([0],axis=-1).flat[:]*(c['x'].take([1],axis=-1).flat[:]-1.0)

        if c.has_key(('velocity',0)):
            c[('velocity',0)][:]=c[('dH',0,0)].flat[:]


T1 = 2.*math.pi

"""
modify to make this on unit square
Example 2
phi_t f(x,y,t) phi_x + g(x,y,t) phi_y = 0

f(x,y,t) = sin^2(pi(x-1))sin(2pi(y-1))cos(t/tau pi)
g(x,y,t) = -sin^2(pi(y-1))sin(2pi(x-1))cos(t/tau pi)

tau = 1.5

phi(x,0) = 0       0.3 <= r
           0.3-r   0.1 < r < 0.3
           0.2     r <= 0.1
r = sqrt((x-05)^2 + (y-0.75)^2)

solution is initial condition at t=tau
periodic bcs on [0,2]^2
"""

L2 = [1.0,1.0,1.0]
class Example2IC:
    def __init__(self,center):
        self.xC = center[0]; self.yC = center[1]
    def uOfXT(self,x,t):
        r = sqrt((x[0]-self.xC)**2 + (x[1]-self.yC)**2)
        if 0.15 <= r:
            return 0.0
        if 0.05 < r and r < 0.15:
            return 0.15-r
        if r <= 0.05:
            return 0.1
 
 
   
class Example2HJ(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import unitSquareVortexLevelSetEvaluate
    def __init__(self,period):
        self.period = period
        
        mass={0:{0:'linear'}}
        advection={0:{0:'linear'}}
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
        self.unitSquareVortexLevelSetEvaluate(t,
                                              c['x'],
                                              c[('u',0)],
                                              c[('grad(u)',0)],
                                              c[('m',0)],
                                              c[('dm',0,0)],
                                              c[('f',0)],
                                              c[('df',0,0)],
                                              c[('H',0)],
                                              c[('dH',0,0)])
T2 = 1.5

"""
Example 3
phi_t + (phi_x + phi_y + 1)^2/2 = 0
phi(x,0) = -cos(x+y)
periodic bcs on [0,2pi]

"""
L3 = [2.0*math.pi,2.0*math.pi,1.0]
class Example3IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -math.cos(x[0]+x[1])
#need to compute actual solution
Example3Soln = Example3IC

#still smooth        
#T3 = 0.1
T3 = 1.

class Example3HJ(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import HJBurgersEvaluate
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
        offset = 1.
        self.HJBurgersEvaluate(offset,
                               c[('u',0)],
                               c[('grad(u)',0)],
                               c[('m',0)],
                               c[('dm',0,0)],
                               c[('H',0)],
                               c[('dH',0,0)])
        
        #c[('m',0)][:]=c[('u',0)]
        #c[('dm',0,0)][:]=1.0
        #c[('H',0)].flat[:]=numpy.inner(c[('grad(u)',0)],c[('grad(u)',0)])
        #c[('H',0)] *= 0.5
        #c[('dH',0,0)].flat[:]=c[('grad(u)',0)].flat[:]
        #c[('velocity',0)][:]=1.0
        if c.has_key(('velocity',0)):
            c[('velocity',0)][:]=c[('dH',0,0)].flat[:]


"""
Example 4
phi_t + (phi_x^2 + phi_y^2)^0.5 = 1
on [0,1] x [0,1] / [0.4,0.6] x [0.4,0.6]

phi(x,0) = max(|x-0.5|,|y-0.5|) - 0.1
phi = 0 on [0.4,0.6] x [0.4,0.6]

"""
L4 = [1,1,1.0]
class Example4IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return max(abs(x[0]-0.5),abs(x[1]-0.5)) - 0.1
#need to compute actual solution
Example4Soln = Example4IC

T4 = 1.#10.

Example4HJ = EikonalEquationCoefficients

Example = 4
if Example == 1:
    analyticalSolution = {0:Example1Soln(Example1IC([0.4,0.4]))}
    coefficients = Example1HJ()
    T = T1
    L=L1
elif Example == 2:
    analyticalSolution = {0:Example2IC([0.5,0.75])}#{0:Example2Soln([1.4,1.4])}
    coefficients = Example2HJ(8.)
    T = T2
    L=L2
elif Example == 3:
    analyticalSolution = {0:Example3Soln()}
    coefficients = Example3HJ()
    T = T3
    L = L3
else:
    assert Example == 4
    polyfile = "chengShuSquare"
    analyticalSolution = {0:Example4Soln()}
    coefficients = Example4HJ(rhsval=1.0)
    T = T4
    L = L4

name = 'ls_ChengShu_2d_Ex4_dgp1_Jaffre_nu0_1_3level'

def getDBC4(x,tag):
    if tag == 6:
        #mwf debug
        print "Cheng-Shu getDBC4 x[0]= %s x[1]= %s tag = %s " % (x[0],x[1],tag)
        return lambda x,t: 0.0
coefficients.variableNames=['phi']

#now define the Dirichlet boundary conditions

def getDBC(x,tag):
    pass


if Example == 4:
    dirichletConditions = {0:getDBC4}
else:
    dirichletConditions = {0:getDBC}

#periodic boundary conditions
def getPDBC(x):
    if (x[0] == 0.0 or x[0] == L[0]) and (x[1] == 0.0 or x[1] == L[1]):
        return numpy.array([0.0,0.0,0.0])
    elif x[0] == 0.0 or x[0] == L[0]:
        return numpy.array([0.0,round(x[1],5),0.0])
    elif (x[1] == 0.0 or x[1] == L[1]):# and (0.0 < x[0] and x[0] < 1):
        return numpy.array([round(x[0],5),0.0,0.0])



if Example != 4:
    periodicDirichletConditions = {0:getPDBC}

initialConditions  = {0:analyticalSolution[0]}

def getAFBC(x,tag):
    pass

fluxBoundaryConditions = {0:'outflow'}
advectiveFluxBoundaryConditions =  {}#can't have for dg {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}

## @}
