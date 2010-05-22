from pyadh import *
from pyadh.default_p import *

"""
Nonlinear advection-diffusion-reaction of a shock in 1D.
"""

##\page Tests Test Problems 
# \ref nladr_1c_shock_1d_p.py "Nonlinear advection-diffusion-reaciton of a shock"
#

##\ingroup test
#\file nladr_1c_shock_1d_p.py
#\brief Nonlinear advection-diffusion-reaction of a shock in 1D.
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

a0=1.0e-1
#a0=0.0
#b0=0.0
b0=1.0
A0={0:Numeric.array([[a0]])}
B0={0:Numeric.array([b0])}
C0={0:0.0}
M0={0:1.0}

p0={0:1.0}
q0={0:2.0}
r0={0:1.0}
s0={0:0.0}
t0={0:0.0}


analyticalSolutions = None

coefficients = NonlinearVADR_pqrst(M=M0,A=A0,B=B0,C=C0,
                                   p=p0,q=q0,r=r0,s=s0,t=t0)

#now define the Dirichlet boundary conditions

def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: 1.0
#    if x[0] == 1.0:
#        return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

class ShockIC:
    def __init__(self,dbc):
        self.uLeft = dbc([0.0,0.0,0.0])([0.0,0.0,0.0],0.0)
        self.uRight=0.0
        #self.uRight = dbc([1.0,0.0,0.0])([1.0,0.0,0.0],0.0)
    def uOfXT(self,x,t):
        if x[0] <= 0.0:
            return self.uLeft
        else:
            return self.uRight

initialConditions  = {0:ShockIC(getDBC)}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 1.0
