import math
nd = 2

polyfile = "hull"

DT = 1.0e-2
T = 10.0
nDTout = int((T/DT))
nLevels = 3

inflow = -1.0
waterLevel = 19.0
slope = 0.0
epsFact = 3.0
#epsFact = 0.0
#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5

#gravity
g=[0.0,-9.8]
LT = 3.25000e+01
LR = 2.0e2

def getDBC_p_susan(x):
    #top (open)
    if x[1] == LT:
        return lambda x,t: 0.0

def getDBC_u_susan(x):
#     #botom and sides (no flow)
#     if (x[0] == 0.0 or
#         x[0] == LR or
#         x[1] == 0.0):
#         return lambda x,t: 0.0
    #sides 
    if (x[0] == LR):
        return lambda x,t: inflow
def getDBC_v_susan(x):
#     #botom and sides (no flow)
#     if (x[0] == 0.0 or
#         x[0] == LR or
#         x[1] == 0.0):
#         return lambda x,t: 0.0
    #bottom (no flow)
    if x[1] == 0.0:
        return lambda x,t: 0.0
    if x[0] == LR:
        return lambda x,t: 0.0
    #top no flow
    #elif x[1] == LT:
    #    return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_susan,
                       1:getDBC_u_susan,
                       2:getDBC_v_susan}

def getAFBC_p_susan(x):
    #botom and sides (no flow)
    if (x[0] == LR):
        return lambda x,t: inflow
    if (x[0] == 0.0 or
        x[1] == LT):
        pass
    else: #no flow
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_susan}

diffusiveFluxBoundaryConditions = {0:{}}

class PerturbedSurface_p:
    def __init__(self,waterLevel,slope):
        self.waterLevel=waterLevel
        self.slope=slope
    def uOfXT(self,x,t):
        z = self.waterLevel + (x[0] - 0.5*LR)*self.slope
        if x[1] > z:
            return -(1.0-x[1])*rho_1*g[1]
        else:
            return -((1.0-z)*rho_1 +
                    (z-x[1])*rho_0
                    )*g[1]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(waterLevel,slope),
                     1:AtRest(),
                     2:AtRest()}

