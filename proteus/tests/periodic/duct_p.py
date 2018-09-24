from proteus import *
from proteus.default_p import *
from proteus.mprans import RANS2P
"""
Navier-Stokes flow in a periodic duct with square cross-section
"""
nd = 3

L=(4.0, 1.0, 1.0)

LevelModelType = RANS2P.LevelModel

coefficients = RANS2P.Coefficients(epsFact=0.0,
                                   sigma=0.0,
                                   rho_0=998.2,nu_0=1.004e-6,
                                   rho_1=998.2,nu_1=1.004e-6,
                                   g=[1.0e-1, 0.0, 0.0],
                                   nd=nd,
                                   ME_model=0,
                                   VF_model=None,
                                   LS_model=None,
                                   Closure_0_model=None,
                                   Closure_1_model=None,
                                   epsFact_density=0.0,
                                   stokes=False,
                                   useVF=0.0,
                                   useRBLES=0.0,
                                   useMetrics=1.0,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=100.0,
                                   forceStrongDirichlet=False,
                                   turbulenceClosureModel=0,
                                   NONCONSERVATIVE_FORM=1.0)

T = 100.0
nsave=100
dt_init = 1.0e-3
DT = (T-dt_init)/float(nsave-1)
tnList = [0.0,dt_init]+[dt_init+i*DT for i in range(nsave)]

eps=1.0e-8

def getDBC_pressure_duct(x,flag):
    return None

def getDBC_u_duct(x,flag):
    if x[2] <  eps or x[2] > L[2] - eps:#top and bottom: no slip
        return lambda x,t: 0.0

def getDBC_v_duct(x,flag):
    if x[2] <  eps or x[2] > L[2] - eps:#top and bottom: no slip
        return lambda x,t: 0.0
    
def getDBC_w_duct(x,flag):
    if x[2] <  eps or x[2] > L[2] - eps:#top and bottom: no slip
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_duct,
                       1:getDBC_u_duct,
                       2:getDBC_v_duct,
                       3:getDBC_w_duct}

periodic = True

if periodic:
    def getPDBC(x,flag):
        if (x[0] < eps or x[0] > L[0] - eps) and (x[1] < eps or x[1] > L[1] - eps) and (x[2] < eps or x[2] > L[2] - eps):#x,y,z corner
            return numpy.array([0.0,0.0,0.0])
        elif (x[0] < eps or x[0] > L[0] - eps) and (x[2] < eps or x[2] > L[2] - eps):#x-z edge
            return numpy.array([0.0,round(x[1],5),0.0])
        elif (x[0] < eps or x[0] > L[0] - eps) and (x[1] < eps or x[1] > L[1] - eps):#x-y edge
            return numpy.array([0.0,0.0,round(x[2],5)])
        elif (x[1] < eps) or (x[1] > L[1]-eps):#on front or back
            return numpy.array([round(x[0],5),0.0,round(x[2],5)])
        elif (x[0] < eps) or (x[0] > L[0]-eps):#on inflow or outflow (left/right)
            return numpy.array([0.0,round(x[1],5),round(x[2],5)])

    periodicDirichletConditions = {0:getPDBC,
                                   1:getPDBC,
                                   2:getPDBC,
                                   3:getPDBC}

def getAFBC_p_duct(x,flag):
    if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
        return lambda x,t: 0.0
    else:
        return lambda x,t: 0.0

def getAFBC_u_duct(x,flag):
    if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
        return lambda x,t: 0.0
    else:
        return lambda x,t: 0.0

def getAFBC_v_duct(x,flag):
    if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
        return lambda x,t: 0.0
    else:
        return lambda x,t: 0.0

def getAFBC_w_duct(x,flag):
    if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no flow
        return lambda x,t: 0.0
    else:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                    1:getAFBC_u_duct,
                                    2:getAFBC_v_duct,
                                    3:getAFBC_w_duct}

def getDFBC_duct(x,flag):
    if (x[2] < eps) or (x[2] > L[2] - eps):#top and bottom: no sliop
        return None
    else:
        return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_duct},
                                   2:{2:getDFBC_duct},
                                   3:{3:getDFBC_duct}}
