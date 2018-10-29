from __future__ import absolute_import
from builtins import object
from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
import numpy as np
from vos import *

LevelModelType = VOF.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=2

coefficients = MyCoefficients(
    checkMass=checkMass,
    FCT=ct.FCT,
    LUMPED_MASS_MATRIX=ct.LUMPED_MASS_MATRIX, 
    STABILIZATION_TYPE=ct.STABILIZATION_TYPE, 
    ENTROPY_TYPE=ct.ENTROPY_TYPE, 
    cE=ct.cE, cK=ct.cK,
    num_fct_iter=1)
coefficients.variableNames=['u']

##################
# VELOCITY FIELD #
##################
radius=0.4
xc=0.5
yc=0.5
def velx(X,t):
    if ct.problem==0:
        return 1.0
    elif ct.problem==1:
        r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
        return (xc-X[0])/(r+1E-10)
    else:
        r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
        return (1-2*X[0])/(r+1E-10)*np.maximum(0,r-0.1)
        
def vely(X,t):
    if ct.problem==0:
        return 0.0
    elif ct.problem==1:
        r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
        return (yc-X[1])/(r+1E-10)    
    else:
        r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
        return (1-2*X[1])/(r+1E-10)*np.maximum(0,r-0.1)
    
velocityFieldAsFunction={0:velx,
                         1:vely}

#####################
# INITIAL CONDITION #
#####################
class init_cond(object):
    def uOfXT(self,X,t):
        if ct.problem==0:
            if X[0] >= 0.05+t and X[0] <= 0.15+t:
                return 1.0
            else:
                return 0.0
        elif ct.problem==1:
            r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
            if r <= radius:
                return 0.5
            else:
                return 0.
        else:
            r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
            if 0.3 <= r and r <= 0.4:
                return 0.5
            else:
                return 0.0
#
initialConditions  = {0:init_cond()}
analyticalSolution = {0:init_cond()}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    return lambda x,t: 0.0
dirichletConditions = {0:getDBC}

def zeroadv(x,flag):
    None
    
advectiveFluxBoundaryConditions =  {0:zeroadv}
fluxBoundaryConditions = {0:'outFlow'}
diffusiveFluxBoundaryConditions = {0:{}}
