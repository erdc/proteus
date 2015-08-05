from proteus import *
from proteus.default_p import *
from splashcube import *
from proteus.mprans import RANS2P
from decimal import *

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

Uinf = 1.0

'''
def getDBC_p(x,flag):
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12:
        return lambda x,t: 0.0
    
def getDBC_u(x,flag):
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12 or flag==boundaryTags['bottom'] or x[2] <= 1.0e-12 or flag==boundaryTags['right'] or x[0] >= L[0] - 1.0e-12 or flag==boundaryTags['left'] or x[0] <= 1.0e-12 or flag==boundaryTags['front'] or x[1]<=1.0e-12 or flag==boundaryTags['back'] or x[1]>L[1]-1.0e-12 : 
        return lambda x,t: 0.0 

def getDBC_v(x,flag):
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12:
    #if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12 or flag==boundaryTags['bottom'] or x[2] <= 1.0e-12 or flag==boundaryTags['right'] or x[0] >= L[0] - 1.0e-12 or flag==boundaryTags['left'] or x[0] <= 1.0e-12 or flag==boundaryTags['front'] or x[1]<=1.0e-12 or flag==boundaryTags['back'] or x[1]>L[1]-1.0e-12 : 
        return lambda x,t: 1.0#Uinf*x[2]/L[2]

def getDBC_w(x,flag):
    #if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12:
    if flag == boundaryTags['top'] or x[2] >= L[2] - 1.0e-12 or flag==boundaryTags['bottom'] or x[2] <= 1.0e-12 or flag==boundaryTags['right'] or x[0] >= L[0] - 1.0e-12 or flag==boundaryTags['left'] or x[0] <= 1.0e-12 or flag==boundaryTags['front'] or x[1]<=1.0e-12 or flag==boundaryTags['back'] or x[1]>L[1]-1.0e-12 : 
        return lambda x,t: 0.0
'''
def getDBC_p(x,flag):
    if flag == boundaryTags['back']:
    #if flag == boundaryTags['top'] or flag==boundaryTags['bottom'] or flag==boundaryTags['front'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:
    #if flag==boundaryTags['top'] or flag==boundaryTags['bottom'] or flag==boundaryTags['front'] or flag==boundaryTags['left'] or flag==boundaryTags['right'] or flag==boundaryTags['back']:
        return lambda x,t: 0.0
    
def getDBC_u(x,flag):
    if flag == boundaryTags['top'] or flag==boundaryTags['bottom'] or flag==boundaryTags['front'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:
    #if flag==boundaryTags['top'] or flag==boundaryTags['bottom'] or flag==boundaryTags['front'] or flag==boundaryTags['left'] or flag==boundaryTags['right'] or flag==boundaryTags['back']:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    #if flag == boundaryTags['top']:
    if flag==boundaryTags['top'] or flag==boundaryTags['bottom'] or flag==boundaryTags['front'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:# or flag==boundaryTags['back']:
        #print x[2], L[2]
        #return lambda x,t: x[2]*1.0/L[2]
        #return lambda x,t: Decimal(x[2])/Decimal(L[2])
        if x[2]>L[2]:
          return lambda x,t: 1.0
        else:
          return lambda x,t: x[2]/L[2]   
        #return lambda x,t: 1.0
    #elif flag ==boundaryTags['top'] or x[2] > L[2]-1.0e-12:
    #    return lambda x,t: 1.0

def getDBC_w(x,flag):
    #if flag == boundaryTags['top']:
    if flag == boundaryTags['top'] or flag==boundaryTags['bottom'] or flag==boundaryTags['front'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:
    #if flag==boundaryTags['top'] or flag==boundaryTags['bottom'] or flag==boundaryTags['front'] or flag==boundaryTags['left'] or flag==boundaryTags['right'] or flag==boundaryTags['back']:
        return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
#        if flag !=boundaryTags['back']:
          #return lambda x,t: 0.0
        pass

def getAFBC_u(x,flag):
        pass

def getAFBC_v(x,flag):
        return lambda x,t: 1.0
        #pass

def getAFBC_w(x,flag):
        pass

def getDFBC_u(x,flag):
    if flag==boundaryTags['back']:
        return lambda x,t: 0.0
#        pass
    
def getDFBC_v(x,flag):
   if flag==boundaryTags['back']:
        return lambda x,t: 0.0
#        pass

def getDFBC_w(x,flag):
    if flag==boundaryTags['back']:
        return lambda x,t: 0.0
#        pass

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

class Couette:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        v = Uinf*(x[2])/L[2]
        return v 

'''
initialConditions = {0:PerturbedSurface_p(waterLine_z),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
'''
initialConditions = {0:AtRest(),
                     1:AtRest(),
                     #2:AtRest(),
                     2:Couette(),
                     3:AtRest()}

