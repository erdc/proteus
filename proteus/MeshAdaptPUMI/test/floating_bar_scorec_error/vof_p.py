from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from floating_bar import *
from proteus.mprans import VOF
from proteus import Context
ct = Context.get()

LevelModelType = VOF.LevelModel
if useOnlyVF:
    RD_model = None
    LS_model = None
else:
    RD_model = 3
    LS_model = 2

coefficients = VOF.Coefficients(LS_model=int(ct.movingDomain)+LS_model,
                                V_model=int(ct.movingDomain)+0,
                                RD_model=int(ct.movingDomain)+RD_model,
                                ME_model=int(ct.movingDomain)+1,
                                checkMass=True,
                                useMetrics=useMetrics,
                                epsFact=epsFact_vof,
                                sc_uref=vof_sc_uref,
                                sc_beta=vof_sc_beta,
                                movingDomain=ct.movingDomain)

def getDBC_vof(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 1.0
    elif flag in [ct.boundaryTags['left'], ct.boundaryTags['right']]:
        if ct.speed > 0.0:
            return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag != boundaryTags['top']:
        return lambda x,t: 0.0
    elif flag in [ct.boundaryTags['left'], ct.boundaryTags['right']]:
        if ct.speed > 0.0:
            return None

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class VF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)

initialConditions  = {0:VF_IC()}
