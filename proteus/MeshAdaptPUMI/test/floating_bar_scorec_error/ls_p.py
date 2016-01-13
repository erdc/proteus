from proteus import *
from proteus.default_p import *
from floating_bar import *
from proteus.mprans import NCLS
from proteus import Context
ct = Context.get()

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=int(ct.movingDomain)+0,
                                 RD_model=int(ct.movingDomain)+3,
                                 ME_model=int(ct.movingDomain)+2,
                                 checkMass=False,
                                 useMetrics=useMetrics,
                                 epsFact=epsFact_consrv_heaviside,
                                 sc_uref=ls_sc_uref,
                                 sc_beta=ls_sc_beta,
                                 movingDomain=ct.movingDomain)

def getDBC_ls(x,flag):
    if flag in [ct.boundaryTags['left'], ct.boundaryTags['right']]:
        if ct.speed > 0.0:
            return lambda x,t: x[2]-waterLevel

dirichletConditions = {0:getDBC_ls}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PHI_IC:
    def uOfXT(self,x,t):
        return x[2] - waterLevel

initialConditions  = {0:PHI_IC()}
