from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from dambreak import *
from proteus.mprans import VOF3P

LevelModelType = VOF3P.LevelModel

coefficients = VOF3P.Coefficients(forceStrongConditions=VOF_forceStrongConditions,
                                  EDGE_VISCOSITY=VOF_EDGE_VISCOSITY,
                                  ENTROPY_VISCOSITY=VOF_ENTROPY_VISCOSITY,
                                  FCT=VOF_FCT,
                                  uL=0,
                                  uR=1,
                                  cK=cK,
                                  cMax=cMax,
                                  cE=cE,
                                  LS_model=LS_model,
                                  V_model=V_model,
                                  RD_model=RD_model,
                                  ME_model=VOF_model,
                                  VOS_model=VOS_model,
                                  checkMass=False,
                                  useMetrics=useMetrics,
                                  epsFact=epsFact_vof,
                                  sc_uref=vof_sc_uref,
                                  sc_beta=vof_sc_beta,
                                  movingDomain=movingDomain)

def getDBC_vof(x,flag):
    if flag == boundaryTags['top'] and openTop:
        return lambda x,t: 1.0

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag != boundaryTags['top'] or not openTop:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))

initialConditions  = {0:VOF_IC()}
