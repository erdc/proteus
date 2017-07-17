from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
nd=2

T=5.
L=(10.0,1.0)
g = 9.81
hl=0.005
xc=5
domain = RectangularDomain(L=L)

#This is relevant just when use_second_order_NonFlatB_with_EV_stabilization=True
cE=1
LUMPED_MASS_MATRIX=0
mannings=0.0

use_SUPG=1

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    return 0.*X[0]

##############################
##### INITIAL CONDITIONS #####
##############################
class dam_break_problem_starting_at_t0:
    def __init__(self,hl=1.,xc=0.5):
        self.hl=hl
        self.xc=xc
    def uOfXT(self,X,t):
        import math
        x = X[0]
        if (x <= self.xc):
            h = self.hl
        else:
            h = 0.001
        return h

class dam_break_problem_starting_at_t1:
    def __init__(self,hl=1.,xc=0.5,g=1):
        self.hl=hl
        self.xc=xc
        self.g=g
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = self.xc-math.sqrt(self.g*self.hl)
        xB = self.xc+2*math.sqrt(self.g*self.hl)
        if (0 <= x and x <= xA):
            return hl
        elif (xA < x <= xB):
            return 4/9./self.g*(math.sqrt(g*self.hl)-(x-self.xc)/2.)**2
        else: 
            return 0.

class velX_starting_at_t1:
    def __init__(self,hl=1.,xc=0.5,g=1):
        self.hl=hl
        self.xc=xc
        self.g=g
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = self.xc-math.sqrt(self.g*self.hl)
        xB = self.xc+2*math.sqrt(self.g*self.hl)
        if (0 <= x and x <= xA):
            return 0.
        elif (xA < x <= xB):
            return 2/3.*(x-self.xc+math.sqrt(g*self.hl))
        else: 
            return 0.

class momX_starting_at_t1:
    def __init__(self,hl=1.,xc=0.5,g=1):
        self.hl=hl
        self.xc=xc
        self.g=g
    def uOfXT(self,X,t):
        h = dam_break_problem_starting_at_t1(hl=self.hl,xc=self.xc,g=self.g)
        vel = velX_starting_at_t1(hl=self.hl,xc=self.xc,g=self.g)
        return h.uOfXT(X,t)*vel.uOfXT(X,t)

class Zero:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:dam_break_problem_starting_at_t0(hl=hl,xc=xc),
                     1:Zero(),
                     2:Zero()}
#initialConditions = {0:dam_break_problem_starting_at_t1(hl=hl,xc=xc,g=g),
#                     1:momX_starting_at_t1(hl=hl,xc=xc,g=g),
#                     2:Zero()}

##########################
##### EXACT SOLUTION #####
##########################
class exact_h_starting_at_t0:
    def __init__(self,hl,xc,g):
        self.hl = hl
        self.xc = xc
        self.g = g
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = self.xc-t*math.sqrt(self.g*self.hl)
        xB = self.xc+2*t*math.sqrt(self.g*self.hl)
        if (0 <= x and x <= xA):
            return hl
        elif (xA < x <= xB):
            return 4/9./self.g*(math.sqrt(g*self.hl)-(x-self.xc)/2./t)**2
        else: 
            return 0.

class exact_h_starting_at_t1:
    def __init__(self,hl,xc,g):
        self.hl = hl
        self.xc = xc
        self.g = g
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = self.xc-(t+1)*math.sqrt(self.g*self.hl)
        xB = self.xc+2*(t+1)*math.sqrt(self.g*self.hl)
        if (0 <= x and x <= xA):
            return hl
        elif (xA < x <= xB):
            return 4/9./self.g*(math.sqrt(g*self.hl)-(x-self.xc)/2./(t+1))**2
        else: 
            return 0.

class exact_velx:
    def __init__(self,hl,xc,g):
        self.hl = hl
        self.g = g
    def uOfXT(self,X,t):
        return 0.0

#analyticalSolution = {0:exact_h_starting_at_t0(hl=hl,xc=xc,g=g), 
analyticalSolution = {0:exact_h_starting_at_t1(hl=hl,xc=xc,g=g), 
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
coefficients = SW2DCV.Coefficients(g=g,bathymetry=bathymetry,cE=cE,LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,mannings=mannings)

