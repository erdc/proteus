from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from oneD_advection import *
from proteus.mprans import VOF
name=soname+"_vof"

"""
VOF for 1D smooth and non smooth profiles
"""

LevelModelType = VOF.LevelModel
coefficients = MyCoefficients(epsFact=epsFactHeaviside,checkMass=checkMass,useMetrics=useMetrics,ME_model=0,
                              EDGE_VISCOSITY=EDGE_VISCOSITY,
                              ENTROPY_VISCOSITY=ENTROPY_VISCOSITY,
                              POWER_SMOOTHNESS_INDICATOR=POWER_SMOOTHNESS_INDICATOR, 
                              LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,
                              FCT=FCT,
                              cK=cK,cE=cE,cMax=cMax)

def Heaviside(phi):
    if phi > 0:
        return 1.0
    elif phi < 0:
        return 0.0
    else:
        return 0.5

class init_cond:
    def __init__(self,center=[0.5,0.75,0.5],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        import numpy as np
        #SMOOTH
        return np.exp(-(X[0]-0.5)**2/0.005)
        #return np.tanh((X[0]-0.25)/0.05)+1
        
        # NON-SMOOTH
        #if (X[0] >= 0.25-0.25/2 and X[0] <= 0.25+0.25/2):
        #    return 1.
        #else: 
        #    return 0.
    

    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end init_cond

class ExactSolution:
    def uOfXT(self,x,t):
        import numpy as np
        #SMOOTH
        return np.exp(-(x[0]-t*0-0.5)**2/0.005)
        #return np.tanh((x[0]-t-0.25)/0.05)+1

        # NON-SMOOTH
        #if (x[0]-t >= 0.25-0.25/2 and x[0]-t <= 0.25+0.25/2):
        #    return 1.
        #else: 
        #    return 0.

        
analyticalSolution = {0:ExactSolution()}

def getDBC(x,flag):
    #return lambda x,t: 1.0
    pass
dirichletConditions = {0:getDBC}

#Periodic Boundary Conditions
eps=1.0e-8
def getPDBC(x,flag):
    if x[0] <= eps or x[0] >= (L[0]-eps):
        return numpy.array([0.0,x[1],0.0])
periodicDirichletConditions = {0:getPDBC}

initialConditions  = {0:init_cond(center=[0.5,0.75],radius=0.15)}
fluxBoundaryConditions = {0:'outFlow'}

#cek made no flux since v.n = 0 for this v
def getAFBC(x,flag):
   return lambda x,t: 1.0
advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
