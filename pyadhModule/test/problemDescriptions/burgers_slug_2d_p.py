from pyadh import *
from pyadh.default_p import *
"""
2D, burgers equation with fixed velocity, Riemann ic
"""

## \page Tests Test Problems 
# \ref burgers_slug_2d_p.py "Burgers equation Riemann problem"
#
#030508 Somewhat contrived Burgers equation example with slug initial condition
# low level priority

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

uSlug=2.0; uBack = 1.0
coefficients = ViscousBurgersEqn(v=[1.0/sqrt(2.0),1.0/sqrt(2.0)],nu=0.0,nd=2)
#now define the Dirichlet boundary conditions
#make sure consistent with ic
def dummyBC(x):
    pass
dirichletConditions = {0:dummyBC}

class SlugIC:
    def __init__(self,uslug=1.0,uback=0.0):
        self.uslug = uslug
        self.uback = uback
    def uOfXT(self,x,t):
        if x[0] >= 0.4 and x[0] <= 0.6 and x[1] >= 0.4 and x[1] <= 0.6: 
            return self.uslug
        else:
            return self.uback

analyticalSolution = {0:SlugIC(uslug=uSlug,uback=uBack)}
initialConditions  = {0:SlugIC(uslug=uSlug,uback=uBack)}

fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.5
