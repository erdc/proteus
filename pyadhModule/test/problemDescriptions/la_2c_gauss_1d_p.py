from pyadh import *
from pyadh.default_p import *
"""
Linear advection of two Gaussians
"""

## \page Tests Test Problems 
# \ref la_2c_gauss_1d_p.py "Linear advection of two components"
#

##\ingroup test
#\file la_2c_gauss_1d_p.py
#@{
#
# \brief Linear advecction of two componets. The initial
# conditions are given by a guassian.
#
# The linear advection equation is
#
#\f[
# u_j + \nable (u_j \mathbf{v_j}) = 0,\quad  j=0,1
#\f]
#
# \todo finish la_2c_gauss_1d_p.py doc

nd = 1

a0=0.0
a1=0.0
b0=1.0
b1=-1.0
A0_2c={0:numpy.array([[a0]]),1:numpy.array([[a1]])}
B0_2c={0:numpy.array([b0]),1:numpy.array([b1])}
C0_2c={0:0.0,1:0.0}
M0_2c={0:1.0,1:1.0}


coefficients = LinearVADR_ConstantCoefficients(nc=2,M=M0_2c,A=A0_2c,B=B0_2c,C=C0_2c)

#now define the Dirichlet boundary conditions

def getDBC0(x):
    if x[0] == 0.0:
        return lambda x,t:  0.0
def getDBC1(x):
    if x[0] == 1.0:
        return lambda x,t:  0.0


dirichletConditions = {0:getDBC0,1:getDBC1}

def getPDBC0(x):
    if x[0] == 0.0 or x[0] == 1.0:
        return numpy.array([0.0,0.0,0.0])

getPDBC1 = getPDBC0

periodicDirichletConditions = {0:getPDBC0,1:getPDBC1}
class GaussIC:
    def __init__(self,sigma=1.0/16.,xc=0.2,b=1.0):
        self.sigma= sigma
        self.xc   = xc
        self.b    = b
    def uOfXT(self,x,t):
        xct= self.xc + self.b*t
        d2 = (x[0]-xct)**2
        return exp(-0.5*d2/self.sigma**2)
    
analyticalSolution = {0:GaussIC(xc=0.2,b=b0),1:GaussIC(xc=0.8,b=b1)}

initialConditions  = {0:GaussIC(xc=0.2,b=b0),1:GaussIC(xc=0.8,b=b1)}

#fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}
fluxBoundaryConditions = {0:'noFlow',1:'noFlow'}

def getAFBC(x):
    if x[0] == 0.0:
        return lambda x,t: 0.0
    if x[0] == 1.0:
        return lambda x,t: 0.0
def getDFBC(x):
    if x[0] == 0.0:
        return lambda x,t: 0.0
    if x[0] == 1.0:
        return lambda x,t: 0.0
    
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC,1:getAFBC}

diffusiveFluxBoundaryConditions = {0:{0:getAFBC},1:{1:getAFBC}}

T = 10.0
