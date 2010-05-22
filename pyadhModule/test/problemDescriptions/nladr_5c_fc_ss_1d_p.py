from pyadh import *
from pyadh.default_p import *
"""
Nonlinear advection-diffusion-reaction system of five components, fully coupled coefficients, at equilibrium.
"""

##\page Tests Test Problems 
# \ref nladr_5c_fc_ss_1d_p.py "Nonlinear advection-diffusion-reaciton of five components, fully coupled, at steady state"
#

##\ingroup test
#\file nladr_5c_fc_ss_1d_p.py
#\brief Nonlinear advection-diffusion-reaction system of five components, fully coupled coefficients, at equilibrium.
#
#The governing equations are described in the NonlinearVADR_pqrst class. The domain is the unit interval.
#\todo describe initial/boundary conditions in nladr_5c_fc_ss_1d_p.py

nd = 1


a=1.0
b=1.0

A={0:Numeric.array([[a]]),
   1:Numeric.array([[a]]),
   2:Numeric.array([[a]]),
   3:Numeric.array([[a]]),
   4:Numeric.array([[a]])}
B={0:Numeric.array([b]),
   1:Numeric.array([b]),
   2:Numeric.array([b]),
   3:Numeric.array([b]),
   4:Numeric.array([b])}
C={0:0.0,
   1:0.0,
   2:0.0,
   3:0.0,
   4:0.0}
M={0:0.0,
   1:0.0,
   2:0.0,
   3:0.0,
   4:0.0}

p={0:0,
   1:0,
   2:0,
   3:0,
   4:0}
q={0:2.2,
   1:1.0,
   2:3.0,
   3:2.5,
   4:1.5}
r={0:1.0,
   1:1.0,
   2:1.0,
   3:1.0,
   4:1.0}
s={0:2.5,
   1:1.2,
   2:1.0,
   3:1.3,
   4:4.3}
t={0:2.5,
   1:5.5,
   2:1.5,
   3:3.5,
   4:2.25}

analyticalSolutions = None

coefficients = NonlinearVADR_pqrst_full(nc=5,M=M,A=A,B=B,C=C,
                                        p=p,q=q,r=r,s=s,t=t)

#now define the Dirichlet boundary conditions

def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: 1.0
    if x[0] == 1.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC,1:getDBC,2:getDBC,3:getDBC,4:getDBC}

class LinearIC:
    def __init__(self,dbc):
        self.uLeft = dbc([0.0,0.0,0.0])([0.0,0.0,0.0],0.0)
        self.uRight = dbc([1.0,0.0,0.0])([1.0,0.0,0.0],0.0)
    def uOfXT(self,x,t):
        return x[0]*self.uRight + (1.0-x[0])*self.uLeft

initialConditions = {0:LinearIC(getDBC),
                     1:LinearIC(getDBC),
                     2:LinearIC(getDBC),
                     3:LinearIC(getDBC),
                     4:LinearIC(getDBC)}


fluxBoundaryConditions = {0:'noFlow',1:'noFlow',2:'noFlow',3:'noFlow',4:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{},3:{},4:{}}

