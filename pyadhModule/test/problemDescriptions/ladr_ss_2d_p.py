from pyadh import *
from pyadh.default_p import *
"""
Linear advection-diffusion-reaction at equilibrium in 1D.
"""

## \page Tests Test Problems 
# \ref ladr_ss_2d_p.py "Linear advection-diffusion-reaction at steady state"
#

##\ingroup test
#\file ladr_ss_2d_p.py
#
#\brief Linear advection-diffusion-reaction at equilibrium in 1D.
#\todo finish ladr_ss_2d_p.py doc

nd = 2
#L=(100.0,100.0,0.0)
L=(1.0,1.0,1.0)

#a0=10.0
a0=0.01
#a0=0.005
b0=1.0
A0_1c={0:numpy.array([[a0,0],[0,a0]])}
B0_1c={0:numpy.array([b0,b0])}
C0_1c={0:0.0}
M0_1c={0:0.0}

#ans = AnalyticalSolutions.LinearAD_SteadyState(b=B0_1c[0][0],a=A0_1c[0][0,0])
#analyticalSolution = {0:ans}
#initialConditions = {0:AnalyticalSolutions.LinearAD_SteadyState(b=0.0,a=1.0)}
initialConditions = None

coefficients = LinearVADR_ConstantCoefficients(M=M0_1c,A=A0_1c,B=B0_1c,C=C0_1c)

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    if x[0] <= 1.0e-8 or x[1] <= 1.0e-8:
        return lambda x,t: 1.0 
    if x[0] >= 1.0 - 1.0e-8 or x[1] >= 1.0-1.0e-8:
        return lambda x,t: 0.0
# def getDBC(x,flag):
#     if x[0] <= 0.0 and x[1] < L[1]:
#         return lambda x,t: 1.0 - math.exp(-10.0*(L[1]-x[1])/L[1])
#     if x[0] >= L[0]:
#         return lambda x,t: 0.0
#     if x[1] <= 0.0 and x[0] < L[0]:
#         return lambda x,t: 1.0 - math.exp(-10.0*(L[0]-x[0])/L[0])
#     if x[1] >= L[1]:
#         return lambda x,t: 0.0
# def getDBC(x,flag):
#     if x[0] <= 0.0 and x[1] < L[1]:
#         return lambda x,t: ans.uOfX(x)
#     if x[0] >= L[0]:
#         return lambda x,t: ans.uOfX(x)
#     if x[1] <= 0.0 and x[0] < L[0]:
#         return lambda x,t: ans.uOfX(x)
#     if x[1] >= L[1]:
#         return lambda x,t: ans.uOfX(x)

dirichletConditions = {0:getDBC}

fluxBoundaryConditions = {0:'noFlow'}

def getFlux(x,flag):
    pass

advectiveFluxBoundaryConditions =  {0:getFlux}

diffusiveFluxBoundaryConditions = {0:{0:getFlux}}

