from pyadh import *
from pyadh.default_p import *
from beach_erosion_board_waves import *
from pyadh import NCLS
if useNCLS:
    LevelModelType = NCLS.OneLevelNCLS
"""
The non-conservative level set description of the free surface of a sloshing two-phase flow in a closed box.
"""
## 

##\ingroup test
#\file ls_so_sloshbox_2d_p.py
#
#\brief The non-conservative level set description of the free surface of a sloshing two-phase flow in a closed box.
#
#\todo finish ls_so_sloshbox_2d_p.py doc

if applyCorrection:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1,checkMass=checkMass)
elif applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=2,ME_model=1,checkMass=checkMass)
else:
    coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1,checkMass=checkMass)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}

class PerturbedSurface_vof:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        surfaceNormal = [-sin(self.slopeAngle),cos(self.slopeAngle)]
        signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[1] - self.waterLevel)*surfaceNormal[1]
        if signedDistance > 0.0:
            vof = 1.0
        else:
            vof = 0.0
        return vof

if runWaveProblem:
    initialConditions = {0:Flat_phi()}
else:
    if VOF:
        initialConditions  = {0:GaussianWaterColumn_vof()}
    else:
        initialConditions  = {0:GaussianWaterColumn_phi()}
    
#if nonconsrv:
fluxBoundaryConditions = {0:'outFlow'}
#else:
#    fluxBoundaryConditions = {0:'noFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
