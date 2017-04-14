from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
nd=2

T=1.0
L=(2.0,2.0)
g = 1.0

cE=2
#Parameters for init condition
wl=0.5
island=1.25
xl=0.75
xr=1.25
yl=1.0-0.25
yr=1.0+0.25

shock=True
domain = RectangularDomain(L=L)
bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    x = X[0]
    y = X[1]
    return island*(xl <= x)*(x <= xr)*(yl <= y)*(y <= yr)
    #return 0.*X[0]

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height:
    def __init__(self,wl=0.1,xl=0.25,xr=0.75,yl=0.25,yr=0.75,island=0.05):
        self.wl=wl
    def uOfXT(self,X,t):
        x = X[0]
        y = X[1]
        
        isl = island*(xl <= x)*(x <= xr)*(yl <= y)*(y <= yr)
        
        flat_surface=False

        if (flat_surface):
            if (xl <= x and x <= xr and yl <= y and y <= yr):
                h=wl-island
            else:
                h=wl
        elif(True):
            if (x < xl/2.):
                h = 2*wl
            elif (xl <= x and x <= xr and yl <= y and y <= yr):
                if (wl>=island):
                    h = wl-island
                else:
                    h = 0.
            else:
                h=wl
        else:
            if (x < xl):
                h=wl
            else:
                h=0.
        return h

class Zero:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:water_height(wl,xl,xr,yl,yr,island),
                     1:Zero(),
                     2:Zero()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
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

#########################################
##### CREATE MODEL AND COEFFICIENTS #####
#########################################
bathymetry={0:bathymetry_function}
LevelModelType = SW2DCV.LevelModel
coefficients = SW2DCV.Coefficients(g=g,bathymetry=bathymetry,cE=cE)
