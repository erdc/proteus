from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
nd=2

#Problem parameters
L=(10.0,1.0)
C=0.01
q0=2
mannings=0.02

T=8.0
g = 9.81
domain = RectangularDomain(L=L)
 
#This is relevant just when use_second_order_NonFlatB_with_EV_stabilization=True
cE=1
LUMPED_MASS_MATRIX=1
USE_EV_BASED_ON_GALERKIN=0

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    return -C*X[0]+C*L[0]

##############################
##### INITIAL CONDITIONS #####
##############################
class h_init_condition:
    def uOfXT(self,x,t):
        n2 = mannings**2
        return (n2*q0**2/C)**(3./10)

class momX_init_condition:
    def uOfXT(self,x,t):
        return q0

class momY_init_condition:
    def uOfXT(self,x,t):
        return 0.

class Zero:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:h_init_condition(),
                     1:momX_init_condition(),
                     2:momY_init_condition()}

##########################
##### EXACT SOLUTION #####
##########################
class exact_h:
    def uOfXT(self,X,t):
        n2 = mannings**2
        return (n2*q0**2/C)**(3./10)

class exact_hu:
    def uOfXT(self,X,t):
        return q0

class exact_hv:
    def uOfXT(self,X,t):
        return 0.

#analyticalSolution = {0:exact_h(),
#                      1:exact_hu(),
#                      2:exact_hv()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def getDBC_h(x,flag):
    n2 = mannings**2
    if (x[0] in [0.0]) or flag in [bt['left']]:
        return lambda x,t: (n2*q0**2/C)**(3./10)
    else:
        return None

#note, these are the same for hu and hv so we can cheat and use  this p-file for SW2DCV and SW2D
def getDBC_u(x,flag):
    if (x[0] in [0.0]) or flag in [bt['left']]:
        return lambda x,t: q0
    else:
        return None

def getDBC_v(x,flag):
    return lambda x,t: 0.0
   #if x[1] in [0.0,L[1]] or flag in [bt['front'],bt['back']]:
   #return lambda x,t: 0.0
   #else:
   #return None

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
coefficients = SW2DCV.Coefficients(g=g,
                                   bathymetry=bathymetry,
                                   cE=cE,
                                   LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,
                                   USE_EV_BASED_ON_GALERKIN=USE_EV_BASED_ON_GALERKIN, 
                                   mannings=mannings)

