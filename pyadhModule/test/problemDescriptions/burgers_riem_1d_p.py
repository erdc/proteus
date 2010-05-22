from pyadh import *
from pyadh.default_p import *
"""
1D, burgers equation Riemann ic
"""

## \page Tests Test Problems 
# \ref burgers_riem_1d_p.py "Burgers equation Riemann problem"
#

#030508 1d scalar nonlinear hyperbolic pde with convex flux
# classical example but low level test problem. 

##\ingroup test
#\file burgers_riem_1d_p.py
#@{
#
# \brief Burgers equation Riemann problem.
#
# The equation is
#
#\f[
# u_j + \pd{u^2/2}{x} = 0
#\f]
#
# \todo finish burgers_riem_1d_p.py doc

nd = 1

coefficients = ViscousBurgersEqn(v=[1.0],nu=0.000001,nd=1)

uL = 2.0; uR = 1.0
def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: uL

dirichletConditions = {0:getDBC}

class RiemIC:
    def __init__(self,dbc):
        self.uLeft = dbc([0.0,0.0,0.0])([0.0,0.0,0.0],0.0)
        self.uRight= uR
    def uOfXT(self,x,t):
        if x[0] <= 0.5: #0.0
            return self.uLeft
        else:
            return self.uRight
    
analyticalSolution = {0:RiemIC(getDBC)}
initialConditions  = {0:RiemIC(getDBC)}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.5
