"""
Incompressible Navier-Stokes flow around a cylinder in 2D.
"""
from proteus import *
from proteus.default_p import *
import sys
from cylinder2d import *
from proteus.mprans import RANS2P

bcsTimeDependent = True
LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   rho_0=rho,
                                   nu_0=nu,
                                   rho_1=rho,
                                   nu_1=nu,
                                   g=g,
                                   nd=nd,
                                   LS_model=None,
                                   epsFact_density=epsFact_density,
                                   stokes=False,#useStokes,
                                   useRBLES=0.0,
                                   useMetrics=0.0,
                                   forceStrongDirichlet=True)


from math import cos,pi
def velRamp(t):
    return 0.41**-2*sin(pi*t/8.0)
    # if t < 2.0:
    #     return 0.5*(1.0-cos(0.5*pi*t))
    # else:
    #     return 1.0

def getDBC_p(x,flag):
    if flag in [boundaryTags['right']]: #,boundaryTags['left'],boundaryTags['front'], boundaryTags['back']]:
        return lambda x,t: 0.0
    
def getDBC_u(x,flag):
    if flag in[ boundaryTags['left']]: #,boundaryTags['right']]:
        return lambda x,t: velRamp(t)*6.0*x[1]*(fl_H-x[1])#*1.5*Um*x[1]*(fl_H-x[1])/(fl_H/2.0)**2
    elif flag in [boundaryTags['obstacle'],boundaryTags['front'],boundaryTags['back']]:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag in[ boundaryTags['left']]: #,boundaryTags['right']]:
        return lambda x,t: 0.0
    elif flag in [boundaryTags['obstacle'],boundaryTags['front'],boundaryTags['back']]:
        return lambda x,t: 0.0
    


dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}

def getAFBC_p(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: -velRamp(t)*6.0*x[1]*(fl_H-x[1]) #velRamp(t)*1.5*Um*x[1]*(fl_H-x[1])/(fl_H/2.0)**2
    elif flag == boundaryTags['right']:
        return None #lambda x,t: velRamp(t)*6.0*x[1]*(fl_H-x[1])
    else:
        return lambda x,t: 0.0
                          

def getAFBC_u(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['obstacle'],boundaryTags['front'],boundaryTags['back'],boundaryTags['right']]:
        return None
    else:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
      if flag in [boundaryTags['left'],boundaryTags['obstacle'],boundaryTags['front'],boundaryTags['back'],boundaryTags['right']]:
          return None
      else:
          return lambda x,t: 0.0

def getDFBC_u(x,flag):
  if flag in [boundaryTags['left'],boundaryTags['obstacle'],boundaryTags['front'],boundaryTags['back']]: #,boundaryTags['right']]:
      return None
  else:
      return lambda x,t: 0.0

def getDFBC_v(x,flag):
  if flag in [boundaryTags['left'],boundaryTags['obstacle'],boundaryTags['front'],boundaryTags['back']]: #,boundaryTags['right']]:
      return None
  else:
      return lambda x,t: 0.0





advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0



initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}
