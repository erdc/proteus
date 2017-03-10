from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
nd=2

T=2.0
L=(2.0,2.0)
g = 1.0
nu= 1.0e-8#1.0e-3#0.1
H0=0.0E-14
HH=2.0
shock=True
domain = RectangularDomain(L=L)
bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

useCV=True
if useCV:
    LevelModelType = SW2DCV.LevelModel
    coefficients = SW2DCV.Coefficients(nu=nu,g=g)
else:
    LevelModelType = SW2D.LevelModel
    coefficients = SW2D.Coefficients(nu=nu,g=g)

class HumpIC:
    def __init__(self,Lx,Ly,r,H,H0,shock=False):
        self.shock=shock
        self.xc = Lx/2.0
        self.yc = Ly/2.0
        self.r=r
        self.H=H
        self.H0=H0
    def uOfXT(self,X,t):
        import math
        x = X[0]
        y = X[1]
        if math.sqrt((x - self.xc)**2 + (y - self.yc)**2) < self.r:
            if self.shock:
                h = self.H
            else:
                h = (self.H*
                     (1.0+cos(math.pi*(x-self.xc)/self.r))*
                     (1.0+cos(math.pi*(y-self.yc)/self.r))
                     +self.H0)
        else:
            h = self.H0
        return h

class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

class ConstIC:
    def uOfXT(self,x,t):
        return 1.0

initialConditions = {0:HumpIC(L[0],L[1],0.25*L[0],HH,H0),
                     1:ZeroIC(),
                     #1:ConstIC(),
                     2:ZeroIC()}

def getDBC_h(x,flag):
    return None

#note, these are the same for hu and hv so we can cheat and use  this p-file for SW2DCV and SW2D
def getDBC_u(x,flag):
   if (x[0] in [0.0,L[0]]) or flag in [bt['left'],bt['right']]:
       return lambda x,t: 0.0
   else:
       return None

def getDBC_v(x,flag):
   if x[1] in [0.0,L[1]] or flag in [bt['front'],bt['back']]:
       return lambda x,t: 0.0
   else:
       return None

dirichletConditions = {0:getDBC_h,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None
def getAFBC_v(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

def getDFBC_u(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None

def getDFBC_v(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}




