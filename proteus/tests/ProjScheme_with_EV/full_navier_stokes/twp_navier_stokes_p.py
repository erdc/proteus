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
                                    VOS_model=VOS_model,
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
                                    PSTAB=1.0, 
                                    cE=cE,
                                    cMax=cMax)

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC_u(x,flag):
    return None

def getDBC_v(x,flag):
    return None

def getAFBC_u(x,flag):
    return lambda x,t: 0.

def getAFBC_v(x,flag):
    return lambda x,t: 0.

def getDFBC_u(x,flag):
    # Set grad(u).n
    if (flag==1): # left boundary 
        return lambda x,t: -2*nu*np.cos(x[0])*np.sin(x[1]+t)*(-1.)
    elif (flag==2): # right boundary 
        return lambda x,t: -2*nu*np.cos(x[0])*np.sin(x[1]+t)*(1.)
    elif (flag==3): # bottom boundary 
        return lambda x,t: 0. 
    else: # top boundary 
        return lambda x,t: 0. 

def getDFBC_v(x,flag):
    if (flag==1): # left boundary 
        return lambda x,t: 0.
    elif (flag==2): # right boundary 
        return lambda x,t: 0.
    elif (flag==3): # bottom boundary 
        return lambda x,t: 2*nu*np.cos(x[0])*np.sin(x[1]+t)*(-1.)
    else: # top boundary 
        return lambda x,t: 2*nu*np.cos(x[0])*np.sin(x[1]+t)*(1.)

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}
advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v}
diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v}}

######################
# INITIAL CONDITIONS #
######################
class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

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
                     1:vely_at_t0()}

###############
# FORCE TERMS #
###############
def forcex(X,t):
    x = X[0]
    y = X[1]
    return (np.sin(x)*np.cos(y+t) # Time derivative
            + np.sin(x)*np.cos(x) # Non-linearity
            + 2*nu*np.sin(x)*np.sin(y+t)) # Diffusion

def forcey(X,t):
    x = X[0]
    y = X[1]
    return (-np.cos(x)*np.sin(y+t) # Time derivative
            - np.sin(y+t)*np.cos(y+t) #Non-linearity
            + 2*nu*np.cos(x)*np.cos(y+t)) # Diffusion

forceTerms = {0:forcex,
              1:forcey}

##################
# EXACT SOLUTION #
##################
class velx:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return np.sin(x[0])*np.sin(x[1]+t)
    def duOfXT(self,x,t):
        return [np.cos(x[0])*np.sin(x[1]+t),
                np.sin(x[0])*np.cos(x[1]+t)]

class vely:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return np.cos(x[0])*np.cos(x[1]+t)
    def duOfXT(self,x,t):
        return [-np.sin(x[0])*np.cos(x[1]+t),
                -np.cos(x[0])*np.sin(x[1]+t)]

analyticalSolution = {0:velx(), 
                      1:vely()}

