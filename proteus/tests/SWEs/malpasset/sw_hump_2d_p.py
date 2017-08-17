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
    ("refinement",4,"Level of refinement"),
    ("reflecting_BCs",1,"Use reflecting BCs")
])

nd=2

T=opts.T
nDTout=opts.nDTout

g = 9.81
# PARAMETERS #
mannings=0.033 #0.025 #0.025#m^{-1/3} after Ying etal 2009
cE=5
LUMPED_MASS_MATRIX=1

domain = None#RectangularDomain(L=L)
meshfile = "mal_50sec"
#
#set bathmetry to zero for debugging
#
#domainTmp = Domain.Mesh2DMDomain(meshfile)
#meshTmp = TriangularMesh()
#meshTmp.generateFrom2DMFile(domainTmp.meshfile)
#meshTmp.nodeArray[:,2]=0.0
#meshTmp.writeMeshADH("mal_50sec_flat")
#meshfile = "mal_50sec_flat"
#
#end debug
#

######################
##### BATHYMETRY #####
######################
#
#using mesh z coord
#
##############################
##### INITIAL CONDITIONS #####
##############################

class water_height_at_t0:
    """set the water level to 100m behind the dam and dry elsewhere"""
    def uOfXT(self,X,t):
        x = X[0]
        #LINE 1
        x1 = 4701.18
        y1 = 4143.41
        x2 = 4655.5
        y2 = 4392.1
        m = (y2-y1)/(x2-x1)
        dam1 = m*(x-x1)+y1

        #LINE 2
        x1 = 4655.5
        y1 = 4392.1
        x2 = 4000.0
        y2 = 5500.0
        m = (y2-y1)/(x2-x1)
        dam2 = m*(x-x1)+y1        

        if (X[1] <= dam1 and X[1] <= dam2): 
            return np.maximum(100.0-X[2],0.)
        else:
            return 0.

class Zero:
    """still water conditions"""
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:water_height_at_t0(),
                     1:Zero(),
                     2:Zero()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def getDBC_h(x,flag):
    return None

def getDBC_u(x,flag):
    return None

def getDBC_v(x,flag):
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
    return lambda x,t: 0.0

def getAFBC_v(x,flag):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

def getDFBC_u(x,flag):
    return lambda x,t: 0.0

def getDFBC_v(x,flag):
    return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

#cek todo, add gauges

#########################################
##### CREATE MODEL AND COEFFICIENTS #####
#########################################
LevelModelType = SW2DCV.LevelModel
coefficients = SW2DCV.Coefficients(g=g,
                                   bathymetry=None,
                                   cE=cE,
                                   LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,
                                   mannings=mannings)
