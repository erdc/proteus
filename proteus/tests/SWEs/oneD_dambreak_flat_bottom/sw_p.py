from __future__ import division
from past.utils import old_div
from builtins import object
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain

# ************************************** #
# ********** PROBLEM SPECIFIC ********** #
# ************************************** #
refinement=7
starting_at_t1=True
# domain and parameters
L=(10.0,1.0)
g = 9.81
hl=0.005
xc=5

# For time and outputting
T=5.
nDTout=100

AUTOMATED_TEST=True
if AUTOMATED_TEST:
    refinement=3
    T=5
    nDTout=1
    
######################
##### PARAMETERS #####
######################
# PHYSICAL PARAMETERS #
g = 9.81
LINEAR_FRICTION = 0
mannings = 0.

# NUMERICAL PARAMETERS #
cE = 1
LUMPED_MASS_MATRIX = 0
SSPOrder = 3
runCFL=0.25
useSuperlu = True
triangleFlag = 1
reflecting_BCs = False

##################
##### DOMAIN #####
##################
nd = 2
domain = RectangularDomain(L=L)
bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

################
##### MESH #####
################
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),10)+1
nnz=1
he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)
domain.MeshOptions.triangleOptions = triangleOptions

######################
##### BATHYMETRY #####
######################
def bathymetry(X):
    return 0.*X[0]

##############################
##### INITIAL CONDITIONS #####
##############################
class dam_break_problem_starting_at_t0():
    def uOfXT(self,X,t):
        import math
        x = X[0]
        if (x <= xc):
            h = hl
        else:
            h = 0.000
        return h

class dam_break_problem_starting_at_t1():
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = xc-math.sqrt(g*hl)
        xB = xc+2*math.sqrt(g*hl)
        if (0 <= x and x <= xA):
            return hl
        elif (xA < x <= xB):
            return 4/9./g*(math.sqrt(g*hl)-old_div((x-xc),2.))**2
        else: 
            return 0.

class velX_starting_at_t1():
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = xc-math.sqrt(g*hl)
        xB = xc+2*math.sqrt(g*hl)
        if (0 <= x and x <= xA):
            return 0.
        elif (xA < x <= xB):
            return 2/3.*(x-xc+math.sqrt(g*hl))
        else: 
            return 0.

class momX_starting_at_t1():
    def uOfXT(self,X,t):
        h = dam_break_problem_starting_at_t1()
        vel = velX_starting_at_t1()
        return h.uOfXT(X,t)*vel.uOfXT(X,t)

class Zero():
    def uOfXT(self,x,t):
        return 0.0

if starting_at_t1:
    IC_h  = dam_break_problem_starting_at_t1
    IC_hu = momX_starting_at_t1
    IC_hv = Zero
else:
    IC_h  = dam_break_problem_starting_at_t0
    IC_hu = Zero
    IC_hv = Zero

###############################
##### BOUNDARY CONDITIONS #####
###############################
def getDBC_h(x,flag):
    return None

def getDBC_hu(x,flag):
   if (x[0] in [0.0,L[0]]) or flag in [bt['left'],bt['right']]:
       return lambda x,t: 0.0
   else:
       return None

def getDBC_hv(x,flag):
    return lambda x,t: 0.0

##########################
##### EXACT SOLUTION #####
##########################
class exact_h_starting_at_t0(object):
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = xc-t*math.sqrt(g*hl)
        xB = xc+2*t*math.sqrt(g*hl)
        if (0 <= x and x <= xA):
            return hl
        elif (xA < x <= xB):
            return 4/9./g*(math.sqrt(g*hl)-(x-xc)/2./t)**2
        else: 
            return 0.

class exact_h_starting_at_t1(object):
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = xc-(t+1)*math.sqrt(g*hl)
        xB = xc+2*(t+1)*math.sqrt(g*hl)
        if (0 <= x and x <= xA):
            return hl
        elif (xA < x <= xB):
            return 4/9./g*(math.sqrt(g*hl)-(x-xc)/2./(t+1))**2
        else: 
            return 0.

if starting_at_t1:
    analyticalSolution = {0:exact_h_starting_at_t1(),
                          1:Zero(),
                          2:Zero()}
else:
    analyticalSolution = {0:exact_h_starting_at_t0(hl=hl,xc=xc,g=g), 
                          1:Zero(),
                          2:Zero()}
   
# ************************************************ #
# ********** GENERIC PORTION OF _p FILE ********** #
# ************************************************ #

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = SW2DCV.LevelModel
coefficients = SW2DCV.Coefficients(g=g,
                                   bathymetry={0:bathymetry},
                                   cE=cE,
                                   LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,
                                   LINEAR_FRICTION=LINEAR_FRICTION,
                                   mannings=mannings)


# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: IC_h(),
                     1: IC_hu(),
                     2: IC_hv()}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu,
                       2:getDBC_hv}
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}
advectiveFluxBoundaryConditions =  {0: lambda x,flag: None,
                                    1: lambda x,flag: None,
                                    2: lambda x,flag: None }
diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1: lambda x,flag: None},
                                   2:{2: lambda x,flag: None}}
