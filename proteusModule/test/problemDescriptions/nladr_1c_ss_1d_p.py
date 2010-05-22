from pyadh import *
from pyadh.default_p import *
"""
Nonlinear advection-diffusion-reaction at equilibriumin 1D.
"""

##\page Tests Test Problems 
# \ref nladr_1c_ss_1d_p.py "Nonlinear advection-diffusion-reaction at steady state"
#

##\ingroup test
#\file nladr_1c_ss_1d_p.py
#\brief Nonlinear advection-diffusion-reaction at equilibrium in 1D.
#
#The governing equations are described in the NonlinearVADR_pqrst class. The domain is the unit interval,
#and the initial/boundary conditions are given by
#\f{eqnarray*}
#u(x,0) &=& 1 \mbox{ for } x \leq 1/2 \\
#u(x,0) &=& 0 \mbox{ for } x > 1/2 \\
#u(0,t) &=& 1 \\
#u_x(1,t) &=& 0
#\f}

nd = 1

initialConditions = None

a0=1.5e-3
b0=1.0
A0={0:numpy.array([[a0]])}
B0={0:numpy.array([b0])}
C0={0:0.0}
M0={0:0.0}

p0={0:0}
q0={0:1}
r0={0:2}
s0={0:0.0}
t0={0:0.0}

if (a0 != 0.0 and
    b0/a0 < 1.0/2.5e-2 and
    q0[0]==2 and
    r0[0]==1):
    analyticalSolutions = {0:AnalyticalSolutions.NonlinearAD_SteadyState(b=B0[0][0],
                                                                         q=q0[0],
                                                                         a=A0[0][0,0],
                                                                         r=r0[0])}
else:
    analyticalSolutions = None
coefficients = NonlinearVADR_pqrst(M=M0,A=A0,B=B0,C=C0,
                                   p=p0,q=q0,r=r0,s=s0,t=t0)

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
#    if x[0] == 0.0:
#        return lambda x,t: 1.0
    if x[0] == 1.0:
        return lambda x,t: 0.75

dirichletConditions = {0:getDBC}

fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

def setInflow(x,flag):
    if x[0] ==0.0:
        return lambda x,t: -10.0

advectiveFluxBoundaryConditions =  {0:setInflow}
