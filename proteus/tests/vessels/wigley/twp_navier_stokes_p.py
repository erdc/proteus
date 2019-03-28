from proteus import *
from proteus.default_p import *
from wigley import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
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

coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0 = rho_0,
                                   nu_0 = nu_0,
                                   rho_1 = rho_1,
                                   nu_1 = nu_1,
                                   g=g,
                                   nd=nd,
                                   VF_model=1,
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
                                   movingDomain=movingDomain)

def getDBC_p(x,flag):
    if openTop and flag == boundaryTags['top']:
        return outflowPressure#outflow
    elif openSides and (flag == boundaryTags['front'] or flag == boundaryTags['back']):
        return outflowPressure#outflow
    elif flag == boundaryTags['right']:
        return outflowPressure#outflow

def getDBC_u(x,flag):
    if flag == boundaryTags['left']:
        return twpflowVelocity_u#inflow
    elif smoothBottom == False and flag == boundaryTags['bottom']:
        return twpflowVelocity_u#no slip
    elif flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0#no slip
    
def getDBC_v(x,flag):
    if flag == boundaryTags['left']:
        return twpflowVelocity_v#inflow
    elif smoothBottom == False and flag == boundaryTags['bottom']:
        return twpflowVelocity_v#no slip
    elif flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0#no slip

def getDBC_w(x,flag):
    if flag == boundaryTags['left']:
        return twpflowVelocity_w#inflow
    elif smoothBottom == False and flag == boundaryTags['bottom']:
        return twpflowVelocity_w#no slip
    elif flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0#no slip

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if flag == boundaryTags['left']:
        return twpflowFlux#inflow
    elif flag == boundaryTags['right']:
        if openEnd:
            return None#outflow
        else:
            return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        if openTop:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    elif flag == boundaryTags['front'] or flag == boundaryTags['back']:
        if openSides:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    elif flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0#no flow
    else:
        return lambda x,t: 0.0#no flow

def getAFBC_u(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        if openEnd:
            return None#outflow
        else:
            return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        if openTop:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    elif flag == boundaryTags['front'] or flag == boundaryTags['back']:
        if openSides:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    else:
        return lambda x,t: 0.0#no flow

def getAFBC_v(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        if openEnd:
            return None#outflow
        else:
            return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        if openTop:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    elif flag == boundaryTags['front'] or flag == boundaryTags['back']:
        if openSides:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    else:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        return None#outflow
    elif flag == boundaryTags['top']:
        if openTop:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    elif flag == boundaryTags['front'] or flag == boundaryTags['back']:
        if openSides:
            return None#outflow
        else:
            return lambda x,t: 0.0#no flow
    else:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['top']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['front'] or flag == boundaryTags['front']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['obstacle']:
        return None#weak Dirichlet
    elif smoothBottom == False and flag == boundaryTags['bottom']:
        return None#weak Dirichlet
    else:
        return lambda x,t: 0.0#no flow or outflow

def getDFBC_v(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['top']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['front'] or flag == boundaryTags['front']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['obstacle']:
        return None#weak Dirichlet
    elif smoothBottom == False and flag == boundaryTags['bottom']:
        return None#weak Dirichlet
    else:
        return lambda x,t: 0.0#no flow our outflow

def getDFBC_w(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['top']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['front'] or flag == boundaryTags['front']:
        return lambda x,t: 0.0#outflow
    elif flag == boundaryTags['obstacle']:
        return None#weak Dirichlet
    elif smoothBottom == False and flag == boundaryTags['bottom']:
        return None#weak Dirichlet
    else:
        return lambda x,t: 0.0#no flow our outflow

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class P_IC:
    def uOfXT(self,x,t):
        return twpflowPressure_init(x,t)

class U_IC:
    def uOfXT(self,x,t):
        return twpflowVelocity_u_init(x,t)

class V_IC:
    def uOfXT(self,x,t):
        return twpflowVelocity_v_init(x,t)

class W_IC:
    def uOfXT(self,x,t):
        return twpflowVelocity_w_init(x,t)

initialConditions = {0:P_IC(),
                     1:U_IC(),
                     2:V_IC(),
                     3:W_IC()}
