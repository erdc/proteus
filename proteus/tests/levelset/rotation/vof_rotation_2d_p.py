from __future__ import absolute_import
from builtins import object
from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
try:
    from .rotation2D import *
except:
    from rotation2D import *
from proteus.mprans import VOF
name=soname+"_vof"

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

LevelModelType = VOF.LevelModel

##\ingroup test
#\file vof_rotation_2d_p.py
#
# \todo finish vof_rotation_2d_p.py

if applyRedistancing:
    coefficients = VOF.Coefficients(LS_model=0,V_model=0,RD_model=1,ME_model=2,checkMass=checkMass,
                                   epsFact=epsFact_vof,useMetrics=useMetrics)
elif not onlyVOF:
    coefficients = VOF.Coefficients(LS_model=0,V_model=0,RD_model=None,ME_model=1,checkMass=checkMass,
                                    epsFact=epsFact_vof,useMetrics=useMetrics)
else:
    coefficients = VOF.Coefficients(RD_model=None,ME_model=0,checkMass=checkMass,
                                    epsFact=epsFact_vof,useMetrics=useMetrics)


def Heaviside(phi):
    if phi > 0:
        return 1.0
    elif phi < 0:
        return 0.0
    else:
        return 0.5

class Rotation_phi(object):
    def __init__(self,center=[0.5,0.75,0.5],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1];
        dBubble = self.radius - sqrt(dx**2 + dy**2)
        return smoothedHeaviside(epsFactHeaviside*he,dBubble)#Heaviside(dBubble)
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Rotation_phi

class Rotation_phi_cylinder(object):
    def __init__(self,center=[0.5,0.75,0.5],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1];
        dBubble = self.radius - sqrt(dx**2 + dy**2)
        return smoothedHeaviside(epsFactHeaviside*he,dBubble)#Heaviside(dBubble)
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Rotation_phi

analyticalSolutions = None

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

initialConditions  = {0:Rotation_phi(center=[0.0,0.5],radius=0.25)}

fluxBoundaryConditions = {0:'outFlow'}

#cek made no flux since v.n = 0 for this v
def getAFBC(x,flag):
   return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
