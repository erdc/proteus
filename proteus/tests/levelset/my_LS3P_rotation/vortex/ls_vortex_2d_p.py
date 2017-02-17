from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from vortex2D import *

LevelModelType = NCLS3P.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

class init_cond:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        #return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        return smoothedHeaviside(epsFactHeaviside*he,self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))

analyticalSolution = {0:init_cond(L)}

RD_model=None
coefficients = MyCoefficients(epsFact=epsFactHeaviside,checkMass=checkMass,RD_model=RD_model,useMetrics=useMetrics,
                              EDGE_VISCOSITY=EDGE_VISCOSITY, 
                              ENTROPY_VISCOSITY=ENTROPY_VISCOSITY,
                              POWER_SMOOTHNESS_INDICATOR=POWER_SMOOTHNESS_INDICATOR, 
                              LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX, 
                              FCT=FCT)

coefficients.variableNames=['u']
initialConditions  = {0:analyticalSolution[0]}

#now define the Dirichlet boundary conditions
def getDBC(x,flag):
    #return lambda x,t: 1.0
    pass
 
def zeroInflow(x):
    return lambda x,t: 0.0

dirichletConditions = {0:getDBC}
fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x,flag):
    return lambda x,t: 0.0

#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}
diffusiveFluxBoundaryConditions = {0:{}}
