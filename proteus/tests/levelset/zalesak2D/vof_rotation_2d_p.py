from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from rotation2D import *
from proteus.mprans import VOF
import ls_rotation_2d_p

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

class Zalesak2D_vof:
    def __init__(self,L=[1.0,
                         1.0],
                 center=[0.5,
                         0.5],
                 radius=0.45,
                 zalesak=True):
        self.phi=ls_rotation_2d_p.Zalesak2D(L,center,radius,zalesak)
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFactHeaviside*he,self.phi.uOfXT(x,t))

analyticalSolutions = {0:Zalesak2D_vof(L=L,
                                       center=[0.0,
                                            0.5],
                                       radius=0.25,
                                       zalesak=True)}

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolutions[0]}

fluxBoundaryConditions = {0:'outFlow'}

def getAFBC(x,flag):
   return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
