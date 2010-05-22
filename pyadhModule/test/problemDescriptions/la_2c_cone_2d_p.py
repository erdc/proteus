from pyadh import *
from pyadh.default_p import *
from math import *

nd = 2
## \page Tests Test Problems 
# \ref la_2c_cone_2d_p.py "Linear advection of two components in a rotating velocity field"
#

##\ingroup test
#\file la_2c_cone_2d_p.py
#@{
#
# \brief Linear advecction of two componets in a rotating velocity field. The initial
# conditions are given by smooth cone shaped functions
#
# The linear advection equation is
#
#\f[
# u_j + \nable (u_j \mathbf{v_j}) = 0,\quad  j=0,1
#\f]
#
# \todo finish la_2dc_cone_2d_p.py doc
class RotatingCone2D:
    def __init__(self,radius):
        self.radius = radius
    def uOfXT(self,x,t):
        centerX = 0.25*sin(2*pi*t) + 0.5
        centerY = 0.25*cos(2*pi*t) + 0.5
        coneX = x[0] - centerX
        coneY = x[1] - centerY
        if sqrt(coneX**2 + coneY**2) < self.radius:
            return 0.25*(1.0 + cos(pi*coneX/self.radius))*(1.0+cos(pi*coneY/self.radius))
        else:
            return 0.0
        
class RotatingCone2D_2:
    def __init__(self,radius):
        self.radius = radius
    def uOfXT(self,x,t):
        centerX = 0.25*sin(2*pi*t+pi) + 0.5
        centerY = 0.25*cos(2*pi*t+pi) + 0.5
        coneX = x[0] - centerX
        coneY = x[1] - centerY
        if sqrt(coneX**2 + coneY**2) < self.radius:
            return 0.25*(1.0 + cos(pi*coneX/self.radius))*(1.0+cos(pi*coneY/self.radius))
        else:
            return 0.0

analyticalSolution = {0:RotatingCone2D(1.0/8.0),1:RotatingCone2D_2(1.0/8.0)}

class UnitSquareRotation2(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import unitSquareRotationEvaluate
    def __init__(self):
        mass={0:{0:'linear'},1:{1:'linear'}}
        advection={0:{0:'linear'},1:{1:'linear'}}
        diffusion={}
        potential={}
        reaction={}
#         diffusion={0:{0:{0:'constant'}},1:{1:{1:'constant'}}}
#         potential={0:{0:'u'},1:{1:'u'}}
#         reaction={0:{0:'linear'},1:{1:'linear'}}
        hamiltonian={}
        TransportCoefficients.TC_base.__init__(self,
                                               2,
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
        self.unitSquareRotationEvaluate(c['x'],
                                        c[('u',1)],
                                        c[('m',1)],c[('dm',1,1)],
                                        c[('f',1)],c[('df',1,1)])
        #c[('a',1,1)].flat[::4] = 1.0
        #c[('a',1,1)].flat[3::4] = 1.0
        #c[('f',1)]*=-1.0
        #c[('df',1,1)]*=-1.0

coefficients = UnitSquareRotation2()

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    
dirichletConditions = {0:getDBC,1:getDBC}

initialConditions  = {0:analyticalSolution[0],1:analyticalSolution[1]}

fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}

def zeroadv(x):
    return lambda x,t: 0.0
#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}

diffusiveFluxBoundaryConditions = {}

T = 1.0
