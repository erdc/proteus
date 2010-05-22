from pyadh import *
from pyadh.default_p import *
"""
Linear advection-diffusion-reaction at equilibrium in 1D.
"""

## \page Tests Test Problems 
# \ref ladr_ss_1d_p.py "Linear advection-diffusion-reaction at steady state"
#

##\ingroup test
#\file ladr_ss_1d_p.py
#
#\brief Linear advection-diffusion-reaction at equilibrium in 1D.
#\todo finish ladr_ss_1d_p.py doc

nd = 1


a0=1.0
#a0=0.01
a0=0.005
b0=1.0
A0_1c={0:numpy.array([[a0]])}
B0_1c={0:numpy.array([b0])}
C0_1c={0:0.0}
M0_1c={0:0.0}

analyticalSolution = {0:AnalyticalSolutions.LinearAD_SteadyState(b=B0_1c[0][0],a=A0_1c[0][0,0])}
#initialConditions = {0:AnalyticalSolutions.LinearAD_SteadyState(b=0.0,a=1.0)}
initialConditions = None

coefficients = LinearVADR_ConstantCoefficients(M=M0_1c,A=A0_1c,B=B0_1c,C=C0_1c)

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    if x[0] <= 0.0+1.0e-16:
        return lambda x,t: 1.0
    if x[0] >= 1.0-1.0e-16:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T=1.0
