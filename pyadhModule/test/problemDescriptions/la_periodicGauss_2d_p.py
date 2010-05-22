from pyadh import *
from pyadh.default_p import *
from math import *
from periodicGauss import *
"""
2D, Linear advection of a guassian with periodic bcs
"""

##  \page Tests Test Problems 
# \ref la_periodicGauss_2d_p.py "Linear advection of a Gaussian"
# \addtogroup test
#
#  \file la_periodicGauss_2d_p.py
# @{
#  

##\ingroup test
#  \brief Conservative linear advection of a cone in a constant
#  velocity field.
#



class ConstantVelocityGaussian2D:
    def __init__(self,sigma=1./8.,b=[1.0,0.0],xc=0.25,yc=0.5):
        self.sigma = sigma
        self.xc= xc
        self.yc= yc
        self.b = b
    def uOfXT(self,x,t):
        centerX = (self.xc + self.b[0]*t)%1.0
        centerY = (self.yc + self.b[1]*t)%1.0
        d2 = (x[0]-centerX)**2 + (x[1]-centerY)**2
        return exp(-0.5*d2/self.sigma**2)

#
analyticalSolution = {0:ConstantVelocityGaussian2D(sigma=sigma,b=velocity,xc=center[0],
                                                   yc=center[1])}

M = {0:1.0}
A = {0:numpy.zeros((nd,nd),'d')}
B = {0:velocity}
C = {0:0.0}

if useHJ:
    coefficients = ConstantVelocityLevelSet(b=velocity)
else:
    coefficients = LinearVADR_ConstantCoefficients(nc=1,M=M,A=A,B=B,C=C)

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x,tag):
    return None
#and periodic boundary conditions
def getPDBC(x,tag):
    if (x[0] == 0.0 or x[0] == 1.0) and (x[1] == 0.0 or x[1] == 1.0):
        return numpy.array([0.0,0.0,0.0])
    elif x[0] == 0.0 or x[0] == 1.0:
        return numpy.array([0.0,round(x[1],5),0.0])
    elif (x[1] == 0.0 or x[1] == 1.0):# and (0.0 < x[0] and x[0] < 1):
        return numpy.array([round(x[0],5),0.0,0.0])
periodicDirichletConditions = {0:getPDBC}
parallelPeriodic=True#False    
dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}


advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

## @}
