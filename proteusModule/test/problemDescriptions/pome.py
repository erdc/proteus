nd = 2

#L = (1.0,1.0,1.0)
#L = (1.0e-1,1.0e-1,1.0)
#L = (1.0e-2,1.0e-2,1.0)
L = (1.0e-3,1.0e-3,1.0)
#L = (2.0,2.0,1.0)
#polyfile="pome_haines_cube"
#L = (1.0e-2,3.0*1.0e-2)
waterLevel = 0.25*L[1]
waterLevel_final = 0.25*L[1]
nn=41
nLevels=1
bubble_quad_order = 3
runCFL = 10.0
T= 3.0e-1
DT = 1.0e-4
T=1.0e-2
DT = 1.0e-6
#DT = 1.0e-5
#DT = 1.0e-5
nDTout = int(round(T/DT))

epsFact = 1.5
#epsFact = 0.0
#water
sigma_01 = 72.8e-3
sigma_01 = 0.0
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1=1.500e-5
#
#rho_1 = rho_0
rho_1 = rho_0
nu_1 = nu_0
#gravity
g=[0.0,-9.8]
#g=[0.0,0.0]

bubbleRadius = 0.15*L[0]
bubbleCenter_x = 0.5*L[0]
bubbleCenter_y = 0.25*L[1]
#bubbleRadius = 0.45*L[0]
#bubbleCenter_y = 0.5*L[1]

def getDBC_p_bubble(x):
    import math
    #pass
    #top (no flow)
#    if x[1] == 1.0:
#       if x[0] == 0.0:
#           return lambda x,t: 0.0
#     #capillary tube
    #top open
    if x[1] == L[1]:
        return lambda x,t: 0.0
    #bottom open
    elif x[1] == 0.0:
        return lambda x,t: -(L[1] - waterLevel_final)*rho_1*g[1] - waterLevel_final*rho_0*g[1]

EPS = 1.0e-8
def getDBC_u_bubble(x):
#     pass
    #sides (no flow)
    if (x[0] == 0.0 or
        x[0] == L[0]):
        return lambda x,t: 0.0
    #top and bottom (no slip)
#     if (x[1] == 0.0 or
#         x[1] == L[1]):
#         return lambda x,t: 0.0
#     if x[1] == 0.0:
#         return lambda x,t: 0.0

def getDBC_v_bubble(x):
    pass
    #sides (no slip)
    if (x[0] == 0.0 or
        x[0] == L[0]):
        return lambda x,t: 0.0
    #top and bottom (no flow)
#     if (x[1] == 0.0 or
#         x[1] == L[1]):
#         return lambda x,t: 0.0
    #capilary tube
    #sides no slip
#    if (x[0] == 0.0 or
#        x[0] == L[0]):
#        return lambda x,t: 0.0
# #     #bottom no flow
# #     if x[1] == 0.0:
# #         return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_bubble,
                       1:getDBC_u_bubble,
                       2:getDBC_v_bubble}

def getAFBC_p_bubble(x):
    #botom and sides (no flow)
    if (x[0] == 0.0 or
        x[0] == L[0]):
        return lambda x,t: 0.0
    
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_bubble}

diffusiveFluxBoundaryConditions = {0:{}}

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        import math
        if x[1] > waterLevel:
            return -(L[1]-x[1])*rho_1*g[1]
        else:
            return -(L[1] - waterLevel)*rho_1*g[1] - (waterLevel-x[1])*rho_0*g[1]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Hydrostatic_p(),
                     1:AtRest(),
                     2:AtRest()}

