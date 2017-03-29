from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
nd=2

import math

L=(4.0,4.0)
g = 9.81
# PARAMETERS #
h0=0.1
a=1.0
domain = RectangularDomain(L=L)
eta = 0.5

omega = math.sqrt(2*g*h0)/a
T=3*2*math.pi/omega

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    import numpy as np
    x = X[0]
    y = X[1] 
    r = np.sqrt((x-L[0]/2.)**2+(y-L[0]/2.)**2)
    return -h0*(1-(r/a)**2)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0:
    def __init__(self,Lx=4.,eta=0.5,h0=0.1,a=1.):
        self.Lx=Lx
        self.eta=eta
        self.h0=h0
        self.a=a
    def uOfXT(self,X,t):
        import numpy as np
        x = X[0]
        y = X[1]

        Lx = self.Lx
        a  = self.a
        h0 = self.h0
        eta = self.eta

        #bathymetry
        r = np.sqrt((x-Lx/2.)**2+(y-Lx/2.)**2)
        z = -h0*(1-(r/a)**2)
        
        #water height at t=0
        h = max(0.,eta*h0/a/a*(2*(x-Lx/2.)) - z)
        return h
        #max(0.,eta*h0/a/a*(2*(x-Lx/2.)*np.cos(omega*t)+2*(y-Lx/2.)*np.sin(omega*t)-eta)-z)

class momX_at_t0:
    def uOfXT(self,X,t):
        return 0.0

class momY_at_t0:
    def __init__(self,g=9.81,Lx=4.,eta=0.5,h0=0.1,a=1.0):
        self.g=g
        self.Lx=Lx
        self.eta=eta
        self.h0=h0
        self.a=a
    def uOfXT(self,X,t):
        import numpy as np
        # COMPUTE WATER HEIGHT
        x = X[0]
        y = X[1]

        g = self.g
        Lx = self.Lx
        a  = self.a
        h0 = self.h0
        eta = self.eta
        
        #bathymetry
        r = np.sqrt((x-Lx/2.)**2+(y-Lx/2.)**2)
        z = -h0*(1-(r/a)**2)
        
        #water height at t=0
        h = max(0.,eta*h0/a/a*(2*(x-Lx/2.)) - z)
    
        # COMPUTE MOMENTUM IN Y
        omega = np.sqrt(2*g*h0)/a
        velY = eta*omega
        return h*velY

class Zero:
    def uOfXT(self,x,t):
        return 0.0

class Const:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:water_height_at_t0(Lx=L[0],eta=eta,h0=h0,a=a),
                     1:momX_at_t0(),
                     2:momY_at_t0(g=g,Lx=L[0],eta=eta,h0=h0,a=a)}

##########################
##### EXACT SOLUTION #####
##########################
class water_height:
    def __init__(self,Lx=4.,eta=0.5,h0=0.1,a=1.):
        self.Lx=Lx
        self.eta=eta
        self.h0=h0
        self.a=a
    def uOfXT(self,X,t):
        import numpy as np
        x = X[0]
        y = X[1]

        Lx = self.Lx
        a  = self.a
        h0 = self.h0
        eta = self.eta

        #bathymetry
        r = np.sqrt((x-Lx/2.)**2+(y-Lx/2.)**2)
        z = -h0*(1-(r/a)**2)

        omega = np.sqrt(2*g*h0)/a        

        #water height at t=0
        #h = max(0.,eta*h0/a/a*(2*(x-Lx/2.)) - z)
        h = max(0.,eta*h0/a/a*(2*(x-Lx/2.)*np.cos(omega*t)+2*(y-Lx/2.)*np.sin(omega*t))-z)
        return h

analyticalSolution = {0:water_height(Lx=L[0],eta=eta,h0=h0,a=a),
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
coefficients = SW2DCV.Coefficients(g=g,bathymetry=bathymetry)


