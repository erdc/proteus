from pyadh import *
from pyadh.default_p import *
"""
2D, Buckley Leverett example from Liu 93 paper
"""

## \page Tests Test Problems 
# \ref Buckley_Leverett_Liu_2d_p.py "Buckley Leverett problem with analytical velocity field"
#
#030508 This problem has 2d scalar nonlinear hyperbolic pde with a nonconvex flux, so it
# has some characteristics that aren't tested by a lot of other test problems
# but they aren't a priority right now (low level test problem). 

##\ingroup test
#\file Buckley_Leverett_Liu_2d_p.py
#@{
#
# \brief Buckley Leverett 5 spot injection example with analytical velocity field.
#
# The equation is
#
#\f[
# u_j + \deld{\vec v f(u)} = 0
# f(u) = u^2/(0.2 - 0.4 u^2 + 1.2 u^2)
# v = (-\pd{p}{x},-\pd{p}{y})^T
# p = 0.01 \log((x^2+y^2)^1/2)
#\f]
#
# \todo finish  doc

nd = 2

polyfile = "fiveSpotCut"

coefficients = BuckleyLeverettLiuExample(nu=0.0,nd=2)
name = "Buckley_Leverett_Liu"
#uI = 0.99; u0 = 0.0
uI = 0.0; u0 = 0.99 #this corresponds to actual example
#now define the Dirichlet boundary conditions
def getDBC(x,tag):
    if x[0] <= 0.05 and x[1] <= 0.05:
        return lambda x,t: uI
#make sure consistent with ic
def dummyBC(x,tag):
    pass
dirichletConditions = {0:getDBC}

#cfl calculation chokes if have zero initial condition
class ConstantIC:
    def __init__(self,u0=0.0,uI=1.0):
        self.u0 = u0
        self.uI = uI
    def uOfXT(self,x,t):
        if x[0] <= 0.05 and x[1] <= 0.05:
            return self.uI
        return self.u0

analyticalSolution = {0:ConstantIC(u0=u0,uI=uI)}
initialConditions  = {0:ConstantIC(u0=u0,uI=uI)}

#try to set advective flux at inflow boundary since velocity is perpendicular to normal
def inflow(x,tag):
    if (x[0] <= 0.05 and x[1] <= 0.05):
        return lambda x,t: -0.1

fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 1.5e1
