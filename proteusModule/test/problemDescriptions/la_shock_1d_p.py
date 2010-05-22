from pyadh import *
from pyadh.default_p import *
"""
1D, Linear advection of a step function.
"""

##\page Tests Test Problems 
# \ref la_shock_1d_p.py "Linear advection of a step function"
#

##\ingroup test
#\file la_shock_1d_p.py
#
#\brief Linear advection of a step funciton in 1D.
#
# The linear advection equation is
#
#\f[
# u_j + \nable (u_j \mathbf{v_j}) = 0,\quad  j=0,1
#\f]
#
# \todo finish la_shock_1d_p.py

nd = 1

a0=0.0
b0=1.0
A0_1c={0:numpy.array([[a0]])}
B0_1c={0:numpy.array([b0])}
C0_1c={0:0.0}
M0_1c={0:1.0}

analyticalSolutions = None

coefficients = LinearVADR_ConstantCoefficients(M=M0_1c,A=A0_1c,B=B0_1c,C=C0_1c)

def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t:  1.0
#    if x[0] == 1.0:
#        return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

class ShockIC:
    def __init__(self,dbc):
        self.uLeft = dbc([0.0,0.0,0.0])([0.0,0.0,0.0],0.0)
        self.uRight = 0.0#dbc([1.0,0.0,0.0])([1.0,0.0,0.0],0.0)
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
