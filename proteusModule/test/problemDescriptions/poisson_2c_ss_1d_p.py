from pyadh import *
from pyadh.default_p import *
"""
Poisson's equations for two components (uncoupled) in 1D
"""

##\page Tests Test Problems 
# \ref poisson_2c_ss_1d_p.py "Poisson's equation for two components"
#

##\ingroup test
#\file poisson_2c_ss_1d_p.py
#
#\brief Heterogensou Poisson's equations for two components (uncoupled) in 1D

nd = 1

initialConditions = None

a0=3.5e-2
a1=3.5e-2
b0=0.0
A0={0:numpy.array([[a0]]),1:numpy.array([[a1]])}
B0={0:numpy.array([b0]),1:numpy.array([b0])}
C0={0:0.0,1:0.0}
M0={0:0.0,1:0.0}        

if a0 > 0.0 and (b0/a0 < 1.0/2.5e-2):
    analyticalSolutions = {0:AnalyticalSolutions.LinearAD_SteadyState(b=B0[0][0],a=A0[0][0,0]),
                           1:AnalyticalSolutions.LinearAD_SteadyState(b=B0[1][0],a=A0[1][0,0])}
else:
    analyticalSolutions = None

coefficients = LinearVADR_ConstantCoefficients(nc=1,M=M0,A=A0,B=B0,C=C0)
coefficients.variableNames=['u0','u1']

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    if x[0] == 0.0:
        return lambda x,t: 1.0
    if x[0] == 1.0:
        return lambda x,t: 0.0
def getDBC2(x,flag):
    if x[0] == 0.0:
        return lambda x,t: 2.0
    if x[0] == 1.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC,1:getDBC2}

fluxBoundaryConditions = {0:'noFlow',1:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

