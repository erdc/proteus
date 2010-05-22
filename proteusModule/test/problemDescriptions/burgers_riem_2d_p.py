from pyadh import *
from pyadh.default_p import *
"""
2D, burgers equation with fixed velocity, Riemann ic
"""

## \page Tests Test Problems 
# \ref burgers_slug_2d_p.py "Burgers equation Riemann problem"
#

##\ingroup test
#\file burgers_slug_2d_p.py
#@{
#
# \brief Burgers equation Riemann problem.
#
# The equation is
#
#\f[
# u_j + \deld{\vec v u^2/2} = 0
#\f]
#
# \todo finish burgers_slug_2d_1d_p.py doc

nd = 2


coefficients = ViscousBurgersEqn(v=[1.0/sqrt(2.0),-1.0/sqrt(2.0)],nu=0.0,nd=2)

uL = 2.0; uR = 1.0
#now define the Dirichlet boundary conditions
def getDBC(x):
    if x[0] == 0.0 or x[1] == 1.0:
        return lambda x,t: uL
#make sure consistent with ic
def dummyBC(x):
    pass
dirichletConditions = {0:dummyBC}

class RiemIC:
    def __init__(self,uleft=1.0,uright=0.0):
        self.uLeft = uleft
        self.uRight= uright
    def uOfXT(self,x,t):
        if x[1]-x[0] >= 0.0:
           return self.uLeft
        else:
            return self.uRight

analyticalSolution = {0:RiemIC(uleft=uL,uright=uR)}
initialConditions  = {0:RiemIC(uleft=uL,uright=uR)}


fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.25#0.5
