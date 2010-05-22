from pyadh import *
from pyadh.default_p import *
"""
Linear advection-diffusion of a shock in 1D.
"""
## \page Tests Test Problems 
# \ref lad_shock_1d_p.py "Linear adveciton-diffusion of a shock"
#

##\ingroup test
#\file lad_shock_1d_p.py
#
#\brief Linear advection-diffusion of a shock in 1D.
#
# The linear advection equation is
#
#\f[
# u_j + \nable (u_j \mathbf{v_j}) = 0,\quad  j=0,1
#\f]
#
# \todo finish lad_shock_1d_p.py doc
name="LinearAdvectionDiffusion"

nd = 1

a0=1.0e-3#1.5e-4
b0=1.0
A0_1c={0:numpy.array([[a0]])}
B0_1c={0:numpy.array([b0])}
C0_1c={0:0.0}
M0_1c={0:1.0}

analyticalSolutions = None

coefficients = LinearVADR_ConstantCoefficients(M=M0_1c,A=A0_1c,B=B0_1c,C=C0_1c)
coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x):
     if x[0] == 0.0:
          return lambda x,t: 1.0
     #if x[0] == 1.0:
     #     return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

class ShockIC:
    def __init__(self,dbc):
        pass
#        self.uLeft = dbc([0.0,0.0,0.0])([0.0,0.0,0.0],0.0)
#         self.uRight = dbc([1.0,0.0,0.0])([1.0,0.0,0.0],0.0)
    def uOfXT(self,x,t):
        if x[0] <= 0.0:
            return 1.0
            #return self.uLeft
        else:
            return 0.0
#             return self.uRight

initialConditions  = {0:ShockIC(getDBC)}

fluxBoundaryConditions = {0:'outFlow'}

def getAFBC(x):
     pass
#     if x[0] == 0.0:
#         return lambda x,t: -1.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.5
