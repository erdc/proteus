from proteus import *
from proteus.default_p import *
from dambreak import *
from proteus.mprans import RANS2P

LevelModelType = RANS2P.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2
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
                                   epsFact_density=epsFact_density,
                                   stokes=False,
                                   useVF=useVF,
				   useRBLES=useRBLES,
				   useMetrics=useMetrics,
                                   eb_adjoint_sigma=1.0,
                                   forceStrongDirichlet=0,
                                   turbulenceClosureModel=ns_closure)

def getDBC_p(x,flag):
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12:
        return lambda x,t: 0.0
    
def getDBC_u(x,flag):
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12:
        return lambda x,t: 0.0

def getDBC_w(x,flag):
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12:
        return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if flag != boundaryTags['top'] or x[2] < L[2] - 1.0e-12:
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag != boundaryTags['top'] or x[2] < L[2] - 1.0e-12:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag != boundaryTags['top'] or x[2] < L[2] - 1.0e-12:
        return lambda x,t: 0.0

def getAFBC_w(x,flag):
    if flag != boundaryTags['top'] or x[2] < L[2] - 1.0e-12:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    return lambda x,t: 0.0
    
def getDFBC_v(x,flag):
    return lambda x,t: 0.0

def getDFBC_w(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

class PerturbedSurface_p:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(L[2] - self.waterLevel)*rho_1*g[2] - (self.waterLevel - x[2])*rho_0*g[2]
        else:
            return -(L[2] - self.waterLevel)*rho_1*g[2]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(waterLine_z),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
