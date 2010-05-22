from pyadh import *
from pyadh.default_p import *
"""
Linear advection-diffusion-reaction system with two-components at equilibrium in 1D.
"""

## \page Tests Test Problems 
# \ref ladr_2c_lw_ss_1d_p.py "Linear advection-diffusion-reaction, two-component, lower triangular system, at steady state"
#

##\ingroup test
#\file ladr_2c_lw_ss_1d_p.py
#
#\brief Linear advection-diffusion-reaction system with two-components at equilibrium in 1D.
#\todo finish ladr_2c_lw_ss_1d_p.py doc

nd = 1

initialConditions = None

a0=1.0
b0=1.0
A0={0:Numeric.array([[a0]]),1:Numeric.array([[a0]])}
B0={0:Numeric.array([b0]),1:Numeric.array([-1.0*b0])}
C0={0:0.0,1:0.0}
M0={0:0.0,1:0.0}        

if a0 > 0.0 and (b0/a0 < 1.0/2.5e-2):
    analyticalSolutions = {0:AnalyticalSolutions.LinearAD_SteadyState(b=B0[0][0],a=A0[0][0,0]),
                           1:AnalyticalSolutions.LinearAD_SteadyState(b=B0[1][0],a=A0[1][0,0])}
else:
    analyticalSolutions = None

coefficients = LinearVADR_ConstantCoefficients_lower(nc=2,M=M0,A=A0,B=B0,C=C0)

#now define the Dirichlet boundary conditions

def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: 1.0
    if x[0] == 1.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC,1:getDBC}

fluxBoundaryConditions = {0:'noFlow',1:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

