from proteus import *
from proteus.default_p import *
from NS_convergence import *
from proteus.mprans import RANS3PF

LevelModelType = RANS3PF.LevelModel
LS_model = None
Closure_0_model = None
Closure_1_model = None

coefficients = RANS3PF.Coefficients(epsFact=epsFact_viscosity,
                                    sigma=0.0,
                                    rho_0 = rho_0,
                                    nu_0 = nu_0,
                                    rho_1 = rho_1,
                                    nu_1 = nu_1,
                                    g=g,
                                    nd=nd,
                                    ME_model=V_model,
                                    PRESSURE_model=PRESSURE_model,
                                    SED_model=SED_model,
                                    VOF_model=VOF_model,
                                    LS_model=LS_model,
                                    Closure_0_model=Closure_0_model,
                                    Closure_1_model=Closure_1_model,
                                    epsFact_density=epsFact_density,
                                    stokes=False,
                                    useVF=useVF,
                                    useRBLES=useRBLES,
                                    useMetrics=useMetrics,
                                    eb_adjoint_sigma=1.0,
                                    eb_penalty_constant=weak_bc_penalty_constant,
                                    forceStrongDirichlet=ns_forceStrongDirichlet,
                                    turbulenceClosureModel=ns_closure,
                                    movingDomain=movingDomain,
                                    dragAlpha=dragAlpha,
                                    PSTAB=1.0)

def getDBC_u(x,flag):
    return lambda x,t: np.sin(x[0])*np.sin(x[1]+t)

def getDBC_v(x,flag):
    return lambda x,t: np.cos(x[0])*np.cos(x[1]+t)

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

def getAFBC_u(x,flag):
    return None

def getAFBC_v(x,flag):
    return None

def getDFBC_u(x,flag):
    return lambda x,t: 0.0

def getDFBC_v(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v}}

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

# EXACT SOLUTION #
class velx_at_t0:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return np.sin(x[0])*np.sin(x[1])

class vely_at_t0:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return np.cos(x[0])*np.cos(x[1])

initialConditions = {0:velx_at_t0(),
                     1:velx_at_t0()}
