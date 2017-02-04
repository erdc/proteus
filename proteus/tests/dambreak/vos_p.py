from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from dambreak import *
from proteus.mprans import VOS3P

LevelModelType = VOS3P.LevelModel

coefficients = VOS3P.Coefficients(LS_model=None,
                                  V_model=V_model,
                                  RD_model=None,
                                  ME_model=0,
                                  checkMass=False,
                                  useMetrics=useMetrics,
                                  epsFact=epsFact_vos,
                                  sc_uref=vos_sc_uref,
                                  sc_beta=vos_sc_beta,
                                  movingDomain=movingDomain)

def getDBC_fos(x,flag):
   pass

dirichletConditions = {0:getDBC_fos}

def getAFBC_fos(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_fos}
diffusiveFluxBoundaryConditions = {0:{}}

class Bed:
    def uOfXT(self,x,t):
        return 0.001 + 0.35*smoothedHeaviside(epsFact_consrv_heaviside*he*3,0.25 - x[1])

initialConditions  = {0:Bed()}
