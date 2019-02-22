from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
import math

opts=Context.Options([
    ("T", 3000, "Length of simulation in seconds"),
    ("nDTout", 3000, "number of time steps to archive"),
    ("reflecting_BCs",1,"Use reflecting BCs")
])

T=opts.T
nDTout=opts.nDTout

######################
##### PARAMETERS #####
######################
# PHYSICAL PARAMETERS #
g = 9.81
LINEAR_FRICTION = 0
mannings=0.033 #0.025 #0.025#m^{-1/3} after Ying etal 2009

# NUMERICAL PARAMETERS #
cE = 1
LUMPED_MASS_MATRIX = 0
SSPOrder = 1
runCFL=0.25
useSuperlu = True
triangleFlag = 1
reflecting_BCs = opts.reflecting_BCs

##################
##### DOMAIN #####
##################
nd=2
domain = None
meshfile = "mal_50sec"
he=None

######################
##### BATHYMETRY #####
######################
#
#using mesh z coord
#

##############################
##### INITIAL CONDITIONS #####
##############################
class IC_h(object):
    """set the water level to 100m behind the dam and dry elsewhere"""
    def uOfXT(self,X,t):
        x = X[0]
        #LINE 1
        x1 = 4701.18
        y1 = 4143.41
        x2 = 4655.5
        y2 = 4392.1
        m = old_div((y2-y1),(x2-x1))
        dam1 = m*(x-x1)+y1

        #LINE 2
        x1 = 4655.5
        y1 = 4392.1
        x2 = 4000.0
        y2 = 5500.0
        m = old_div((y2-y1),(x2-x1))
        dam2 = m*(x-x1)+y1        

        if (X[1] <= dam1 and X[1] <= dam2): 
            return np.maximum(100.0-X[2],0.)
        else:
            return 0.

class Zero(object):
    """still water conditions"""
    def uOfXT(self,x,t):
        return 0.0

IC_hu=Zero
IC_hv=Zero

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def getDBC_h(x,flag):
    return None

def getDBC_hu(x,flag):
    return None

def getDBC_hv(x,flag):
    return None

# ************************************************ #
# ********** GENERIC PORTION OF _p FILE ********** #
# ************************************************ #

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = SW2DCV.LevelModel
coefficients = SW2DCV.Coefficients(g=g,
                                   bathymetry=None,
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
