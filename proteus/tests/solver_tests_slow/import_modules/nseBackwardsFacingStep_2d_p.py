from proteus import *
from proteus.default_p import *
from proteus import Domain
import nseBackwardsFacingStep_2d
from TwophaseNavierStokes_ST_LS_SO_VV import TwophaseNavierStokes_ST_LS_SO_VV
from proteus import Context
Context.setFromModule(nseBackwardsFacingStep_2d)
ct=Context.get()

###############################################

# variables set from context
name = ct.name
numeric_scheme = ct.numeric_scheme

###############################################

# space dimension
nd = 2

# Domain Creation
if numeric_scheme == "TH":
    domain = Domain.PlanarStraightLineGraphDomain(
    vertices   =[[-1.0,0.0],[-1.0,1.0],[5.0,1.0],[5.0,-1.0],[0.0,-1.0],[0.0,0.0]],
    vertexFlags=[ 1        , 1        , 2       , 4        , 4        , 4       ],
    segments    =[[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]],
    segmentFlags=[ 1   , 2   , 3   , 4   , 4   , 4   ],
    name="backwardsFacingStep",
    units='cm')
    domain.writePoly('backwardsFacingStep')
elif numeric_scheme == "THQuads":
    #Need to get the quads mesh through a .mat file coming from IFISS
    domain = Domain.MeshQuadDomain_IFISS('grid_data/grid_data_step_4')
else:
    raise Exception, 'Element type not supported yet'

###############################################

eps = 1.0e-8
he = 0.01
#he *=0.5
#he *=0.5
#he *=0.5

vent=False

class uTrue:
    def __init__(self):
        pass
    def uOfX(self,x):
        return 1.
    def uOfXT(self,x,t):
        return 4.0*x[1]*(1-x[1]) 

class uTrue_Re_through_bdy:
    def __init__(self):
        pass
    def uOfX(self,x):
        return 1.
    def uOfXT(self,x,t):
        return (10.0**max(0.0,t-1))*4.0*x[1]*(1-x[1])

class vTrue:
    def __init__(self):
        pass
    def vOFX(self,x):
        return 0.0
    def vOfXT(self,x,t):
        return self.vOFX(x)

def left(x,flag):
    if x[0]==-1.:
        return True
    else:
        return False

def top(x,flag):
    if x[1]==1.:
        return True
    else:
        return False

def right(x,flag):
    if x[0]==5.:
        return True
    else:
        return False

def bottom(x,flag):
    if x[1]==-1:
        return True
    elif (x[0]==0.0 and x[1]<=0.0):
        return True
    elif (x[1]==0.0 and x[0]<=0.0):
        return True
    else:
        return False

def getDBCp(x,flag):
    if right(x,flag):
        return lambda x,t: 0.0
    else:
        pass

def getDBCu(x,flag):
    if left(x,flag):
        return lambda x,t: uTrue().uOfXT(x,t)
    elif top(x,flag) or bottom(x,flag):
        return lambda x,t: 0.0
    else:
        pass

def getDBCu_RE_through_bdy(x,flag):
    if left(x,flag):
        return lambda x,t: utrue_RE_through_bdy().uOfXT(x,t)
    elif top(x,flag) or bottom(x,flag):
        return lambda x,t: 0.0
    else:
        pass

def getDBCv(x,flag):
    if left(x,flag) or top(x,flag) or bottom(x,flag):
        return lambda x,t: 0.0
    else:
        pass

def getAdvFluxBCp(x,flag):
    if left(x,flag):
        return lambda x,t: -uTrue().uOfXT(x,t)
    else:
        pass

def getAdvFluxBCp_RE_through_bdy(x,flag):
    if left(x,flag):
        return lambda x,t: -uTrue_RE_through_bdy().uOfXT(x,t)
    else:
        pass
    
def getAdvFluxBCu(x,flag):
    pass

def getAdvFluxBCv(x,flag):
    pass

def getDiffFluxBCp(x,flag):
    pass

def getDiffFluxBCu(x,flag):
    if right(x,flag):
        return lambda x,t: 0.0
    else:
        pass

def getDiffFluxBCv(x,flag):
    if right(x,flag):
        return lambda x,t: 0.0
    else:
        pass


Diffusivefluxboundaryconditions = {0:getDiffFluxBCp,
                                   1:getDiffFluxBCu,
                                   2:getDiffFluxBCv}

coefficients = TwophaseNavierStokes_ST_LS_SO_VV(epsFact=0.0,
                                                sigma=0.0,
                                                rho_0=1.0,
                                                nu_0=1.0,
                                                rho_1=1.0,
                                                nu_1=1.0,
                                                g=[0.0,0.0],
                                                nd=2,
                                                LS_model=None,
                                                KN_model=None,
                                                epsFact_density=None,
                                                stokes=False);    

coefficients.variableNames=['p','u','v']
