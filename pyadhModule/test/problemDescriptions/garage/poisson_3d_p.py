from pyadh import *
from pyadh.default_p import *
"""
Poisson's equations for two components (uncoupled) in 2D
"""

##\page Tests Test Problems 
# \ref poisson_2d_p.py "Poisson's equation for two couponents"
#

##\ingroup test
#\file poisson_2d_p.py
#
#\brief Poisson's equations for two components (uncoupled) in 2D

nd = 3

initialConditions = None

a0=1.0#3.5e-3
b0=[0.0,0.0,0.0]
A0={0:Numeric.array([[a0,0.0,0.0],[0.0,a0,0.0],[0.0,0.0,a0]])}
B0={0:Numeric.array([b0])}
C0={0:0.0}
M0={0:0.0}

class linearSolution:
    def __init__(self,a0=1.0,ax=0):
        self.a0 = a0
        self.ax = ax
    def uOfX(self,X):
        return 1.0-X[self.ax]
    def uOfXT(self,X,T):
        return 1.0-X[self.ax]

analyticalSolution = {0:linearSolution()}
analyticalSolution = {0:AnalyticalSolutions.PoissonsEquation(K=a0,nd=nd)}

coefficients = LinearVADR_ConstantCoefficients(nc=1,M=M0,A=A0,B=B0,C=C0,rFunc=analyticalSolution)

#now define the Dirichlet boundary conditions

def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: 1.0
    if x[0] == 1.0:
        return lambda x,t: 0.0
    #dir bcs everywhere?
    if x[1] == 0.0 or x[1] == 1.0:
        return lambda x,t: 1.0-x[0]
eps = 1.0e-16
def getDBC(x):
    if x[0] <= 0.0-eps or x[0] >= 1.0+eps:
        return lambda x,t: 0.0
    if x[1] >= 0.0-eps or x[1] >= 1.0+eps:
        return lambda x,t: 0.0
    if x[2] >= 0.0-eps or x[2] >= 1.0+eps:
        return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

