from proteus import *
from proteus.default_p import *
from dambreak import *
from proteus.mprans import RANS3PSed

LevelModelType = RANS3PSed.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if useRANS >= 1:
    Closure_0_model = 5; Closure_1_model=6
    if useOnlyVF:
        Closure_0_model=2; Closure_1_model=3
    if movingDomain:
        Closure_0_model += 1; Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

coefficients = RANS3PSed.Coefficients(epsFact=epsFact_viscosity,
                                      sigma=0.0,
                                      rho_0 = rho_s,
                                      nu_0 = nu_s,
                                      rho_1 = rho_s,
                                      nu_1 = nu_s,
                                      g=g,
                                      nd=nd,
                                      ME_model=SED_model,
                                      PRESSURE_model=PRESSURE_model,
                                      FLUID_model=V_model,
                                      VOS_model=VOS_model,
                                      VF_model=VOF_model,
                                      LS_model=LS_model,
                                      Closure_0_model=Closure_0_model,
                                      Closure_1_model=Closure_1_model,
                                      epsFact_density=epsFact_density,
                                      stokes=False,
                                      useVF=True,
                                      useRBLES=useRBLES,
                                      useMetrics=useMetrics,
                                      eb_adjoint_sigma=1.0,
                                      eb_penalty_constant=weak_bc_penalty_constant,
                                      forceStrongDirichlet=ns_forceStrongDirichlet,
                                      turbulenceClosureModel=ns_closure,
                                      movingDomain=movingDomain,
                                      dragAlpha=dragAlpha)

def getDBC_u(x,flag):
    return None

def getDBC_v(x,flag):
    return None

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

def getAFBC_u(x,flag):
    return lambda x,t: 0.0

def getAFBC_v(x,flag):
    return lambda x,t: 0.0

def getDFBC_u(x,flag):
    return lambda x,t: 0.0

def getDFBC_v(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v}}

class PerturbedSurface_vosStar:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] < 0.25 :
            return 0.3
        else:
            return 0.01

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest()}
