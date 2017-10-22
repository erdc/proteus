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
                                    cMax=cMax,
                                    CORRECT_VELOCITY=CORRECT_VELOCITY)

########################
# BOUNDARY CONDITIONS #
#######################
def getDBC_u(x,flag):
    #None
    pi = np.pi
    #YY: uncomment the next line will produce an error since each edge will be treated
    #as a Dirichlet bc. This ``bug'' is hard to be detected since the result show
    #a small error
    #return lambda x,t: np.sin(x[0])*np.sin(x[1]+t)
    if (flag==1 or flag==2 or flag==3 or flag==4):
        return lambda x,t: np.sin(x[0])*np.sin(x[1]+t)

def getDBC_v(x,flag):
    #None
    pi = np.pi
    #YY: Comment the next line will produce an error since each edge will be treated
    #as a Dirichlet bc. This ``bug'' is hard to be detected since the result show
    #a small error
    #return lambda x,t: np.cos(x[0])*np.cos(x[1]+t)
    if (flag==1 or flag==2 or flag==3 or flag==4):
        return lambda x,t: np.cos(x[0])*np.cos(x[1]+t)

def getAFBC_u(x,flag):
    pass
def getAFBC_v(x,flag):
    pass
def getDFBC_u(x,flag):
    pass
def getDFBC_v(x,flag):
    pass

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}
#advectiveFluxBoundaryConditions =  {}
#diffusiveFluxBoundaryConditions = {}

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
        if manufactured_solution == 1:
            return np.sin(x[0])*np.sin(x[1])
        else:
            return 0.

class vely_at_t0:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if manufactured_solution == 1:
            return np.cos(x[0])*np.cos(x[1])
        else:
            return 0.

initialConditions = {0:velx_at_t0(),
                     1:vely_at_t0()}

#############################
# MATERIAL PARAMETER FIELDS #
#############################
def density(X,t):
    x = X[0]
    y = X[1]
    #return np.sin(x+y)**2+1
    return 1.0

mu_constant=False

def dynamic_viscosity(X,t):
    x = X[0]
    y = X[1]
    return mu*(np.cos(x+y+t)**2+1)

materialParameters = {'density':density,
                      'dynamic_viscosity':dynamic_viscosity}


###############
# FORCE TERMS #
###############
def forcex(X,t):
    x = X[0]
    y = X[1]
    rho = density(X,t)
    return rho*np.sin(x)*np.cos(y+t) # Time derivative


def forcey(X,t):
    x = X[0]
    y = X[1]
    rho = density(X,t)
    return -rho*np.cos(x)*np.sin(y+t) # Time derivative

forceTerms = {0:forcex,
              1:forcey}

##################
# EXACT SOLUTION #
##################
class velx:
    def __init__(self):
        pass
    def uOfXT(self,x,t,flag=None):
        pi = np.pi
        if manufactured_solution == 1:
            return np.sin(x[0])*np.sin(x[1]+t)
        else:
            return np.sin(pi*x[0])*np.cos(pi*x[1])*np.sin(t)
    def duOfXT(self,x,t,flag=None):
        if manufactured_solution == 1:
            return [np.cos(x[0])*np.sin(x[1]+t),
                    np.sin(x[0])*np.cos(x[1]+t)]
        else:
            return [pi*np.cos(pi*x[0])*np.cos(pi*x[1])*np.sin(t),
                    -pi*np.sin(pi*x[0])*np.sin(pi*x[1])*np.sin(t)]

class vely:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        pi = np.pi
        if manufactured_solution == 1:
            return  np.cos(x[0])*np.cos(x[1]+t)
        else:
            return -np.cos(pi*x[0])*np.sin(pi*x[1])*np.sin(t)
    def duOfXT(self,x,t):
        if manufactured_solution == 1:
            return [-np.sin(x[0])*np.cos(x[1]+t),
                    -np.cos(x[0])*np.sin(x[1]+t)]
        else:
            return [pi*np.sin(pi*x[0])*np.sin(pi*x[1])*np.sin(t),
                    -pi*np.cos(pi*x[0])*np.cos(pi*x[1])*np.sin(t)]

analyticalSolution = {0:velx(),
                      1:vely()}

class pressure:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return np.cos(x[0])*np.sin(x[1]+t)
analyticalPressureSolution={0:pressure()}
