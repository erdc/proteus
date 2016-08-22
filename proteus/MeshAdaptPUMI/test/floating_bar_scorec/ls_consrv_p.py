from proteus import *
from proteus.default_p import *
from floating_bar import *
from proteus.mprans import MCorr
from proteus import Context
ct = Context.get()

LevelModelType = MCorr.LevelModel

coefficients = MCorr.Coefficients(LSModel_index=int(ct.movingDomain)+2,
                                  V_model=int(ct.movingDomain)+0,
                                  me_model=int(ct.movingDomain)+4,
                                  VOFModel_index=int(ct.movingDomain)+1,
                                  applyCorrection=applyCorrection,
                                  nd=nd,
                                  checkMass=True,
                                  useMetrics=useMetrics,
                                  epsFactHeaviside=epsFact_consrv_heaviside,
                                  epsFactDirac=epsFact_consrv_dirac,
                                  epsFactDiffusion=epsFact_consrv_diffusion)

class zero_phi:
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return 0.0

initialConditions  = {0:zero_phi()}
