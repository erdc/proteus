from pyadh import *
from pyadh.default_p import *
from math import *
"""
2D, Linear advection of a circular level set function
"""

## \page Tests Test Problems 
# \ref la_rotatingcircle_2d_p.py "Linear advection of a circular level set function"
#

##\ingroup test
#\file la_rotatingcircle_2d_p.py
#
#\brief Linear advection of a circular level set function in a
#rotating velocity field.
# The linear advection equation is
#
#\f[
# u_j + \nable (u_j \mathbf{v_j}) = 0,\quad  j=0,1
#\f]
#
# \todo finish la_rotatingcircle_2d_p.py doc

nd = 2

class RotatingVelocityCircle:
    def __init__(self,radius):
        self.radius = radius
    def uOfXT(self,x,t):
        centerX = 0.25*sin(2*pi*t) + 0.5
        centerY = 0.25*cos(2*pi*t) + 0.5
        return sqrt((x[0]-centerX)**2 + (x[1]-centerY)**2)-self.radius

r0 = 1./8.
analyticalSolution = {0:RotatingVelocityCircle(r0)}

class UnitSquareRotation(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import unitSquareRotationEvaluate
    def __init__(self):
        mass={0:{0:'linear'}}
        advection={0:{0:'linear'}}
        diffusion={0:{0:{0:'constant'}}}
        potential={0:{0:'linear'}}
        reaction={0:{0:'linear'}}
        hamiltonian={}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
    def evaluate(self,t,c):
        self.unitSquareRotationEvaluate(c['x'],
                                        c[('u',0)],
                                        c[('m',0)],c[('dm',0,0)],
                                        c[('f',0)],c[('df',0,0)])

coefficients = UnitSquareRotation()

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
#     if (x[0] == 0.0 or
#         x[0] == 1.0 or
#         x[1] == 0.0 or
#         x[1] == 1.0):
#         return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x):
    return lambda x,t: 0.0
def exactadv(x):
    return lambda x,t: analyticalSolution[0].uOfXT(x,t)
advectiveFluxBoundaryConditions =  {}
#need to check NumericalFlux to make sure this only when inflow 
#advectiveFluxBoundaryConditions =  {0:exactadv}

diffusiveFluxBoundaryConditions = {0:{}}

T = 1.0
