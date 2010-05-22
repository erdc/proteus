from pyadh import *
from pyadh.default_p import *
from beach_erosion_board_waves_3d import *
from pyadh import MCorr
if useMCorr:
    LevelModelType = MCorr.OneLevelMCorr
"""
The non-conservative level set description of a bubble in a two-phase flow
"""

##\ingroup test
#\file ls_bubble_3d_p.py
#
# \todo finish ls_bubble_3d_p.py

coefficients = LevelSetConservation(applyCorrection=applyCorrection,
                                    epsFactHeaviside=epsFact_consrv_heaviside,
                                    epsFactDirac=epsFact_consrv_dirac,
                                    epsFactDiffusion=epsFact_consrv_diffusion,
                                    LSModel_index=1,V_model=0,me_model=4,VOFModel_index=2,nd=nd)

#waterLevel = 0.9*L[1]

class zero_phi:
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return self.uOfX(X)


def getDBC_cnsrv(x,flag):
    pass


dirichletConditions = {0:getDBC_cnsrv}
#bubble rise
initialConditions  = {0:zero_phi()}

fluxBoundaryConditions = {0:'outFlow'}


def getAFBC_cnsrv(x,flag):
    pass

def getDFBC_cnsrv(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_cnsrv}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_cnsrv}}
