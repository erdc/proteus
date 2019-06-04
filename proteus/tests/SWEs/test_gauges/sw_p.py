from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import SW2D,SW2DCV
from proteus.Domain import RectangularDomain, PlanarStraightLineGraphDomain

import numpy as np
import math

opts=Context.Options([
    ("T", 10.0, "Length of simulation in seconds"),
    ("nDTout", 100, "number of time steps to archive"),
    ("refinement",3,"Level of refinement"),
    ("structured",False,"Use structured mesh"),
    ("reflecting_BCs",1,"Use reflecting BCs")
],mutable=True)

AUTOMATED_TEST=True
if AUTOMATED_TEST:
    opts.T=10
    opts.nDTout=1
    opts.refinement=2
#
refinement=opts.refinement
T=opts.T
nDTout=opts.nDTout
L=(75.0,30.0)

######################
##### PARAMETERS #####
######################
# PHYSICAL PARAMETERS #
g = 9.81
LINEAR_FRICTION = 0
mannings = 0.02

# NUMERICAL PARAMETERS #
cE = 1
LUMPED_MASS_MATRIX = 0
SSPOrder = 3
runCFL=0.25
useSuperlu = True
triangleFlag = 1
reflecting_BCs = opts.reflecting_BCs

##################
##### DOMAIN #####
##################
nd=2
domainRect = RectangularDomain(L=L)
if opts.structured:
    domain=domainRect
else:
    domainRect.writePoly("hump")
    domain = PlanarStraightLineGraphDomain("hump")
    domain.boundaryTags = domainRect.boundaryTags
bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
if opts.structured:
    domain.writePoly("tank2d")

################
##### MESH #####
################
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),2)+1
nnz = 1
he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)
domain.MeshOptions.triangleOptions = triangleOptions

######################
##### BATHYMETRY #####
######################
def bathymetry(X):
    import numpy as np
    x = X[0]
    y = X[1]
    bump1 = 1-1./8*np.sqrt((x-30)**2+(y-6)**2)
    bump2 = 1-1./8*np.sqrt((x-30)**2+(y-24)**2)
    bump3 = 3-3./10*np.sqrt((x-47.5)**2+(y-15)**2)
    return np.maximum(np.maximum(np.maximum(0.,bump1),bump2),bump3)

##############################
##### INITIAL CONDITIONS #####
##############################
class IC_h(object):
    def uOfXT(self,X,t):
        x = X[0]
        if (x <= 16):
            eta=1.875
        else:
            eta=0.
        z = bathymetry(X)
        return max(eta - z,0.)

class Zero(object):
    def uOfXT(self,x,t):
        return 0.0

IC_hu = Zero
IC_hv = Zero

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
