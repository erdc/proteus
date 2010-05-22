applyCorrection=True
applyRedistancing=True

nd = 3

L = (1.0e-1,1.0e-1,2.0e-1)

nn=12
nnx=11
nny=11
nnz=21
nLevels=1

useStokes=False
epsFact_density = 3.0
epsFact_viscosity = 3.0
epsFact_redistance = 1.5
epsFact_curvature = 3.0
epsFact_consrv_heaviside = 0.1
epsFact_consrv_dirac = 0.1
epsFact_consrv_diffusion=10.0
bubble_quad_order = 3

runCFL = 10.0
T= 1.0
DT=1.0e-3
nDTout = int(round(T/DT))

#Heaviside/Dirac smoothing factor 
epsFact = 3.0

#water
sigma_01 = 72.8e-3
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1=1.500e-5

#gravity
g=[0.0,0.0,-9.8]

#initial bubble shape
#circle
bubbleRadius = 0.25*L[0]
bubbleCenter_x = 0.5*L[0]
bubbleCenter_y = 0.5*L[1]
bubbleCenter_z = 0.25*L[2]
#elongate by factor along principle axes
ellipticalBubble_x = 1.0
ellipticalBubble_y = 1.0
ellipticalBubble_z = 1.0

#top water surface; can be above or below top of domain
waterLevel = 0.75*L[2]

#Navier-stokes
def getDBC_p_bubble(x):
    #top open
    if x[2] == L[2]:
        return lambda x,t: 0.0
#     #No flow except at one upper corner
#     if x[0] == 0.0 and x[1] == 0.0 and x[2] == L[2]:
#         return lambda x,t: 0.0

def getDBC_u_bubble(x):
    #slip
    pass
#     #sides (no flow and no slip)
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0

def getDBC_v_bubble(x):
    #slip
    pass
#     #sides (no flow and no slip)
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0

def getDBC_w_bubble(x):
    #slip
    pass
    #sides (no flow and no slip)
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0

def getAFBC_p_bubble(x):
    #sides and bottom no flow
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1] or
        x[2] == 0.0):
        return lambda x,t: 0.0
    else:#top open
        pass
#     #sides (no flow)
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1] or
#         x[2] == 0.0 or
#         x[2] == L[2]):
#         return lambda x,t: 0.0

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        import math
        pure_water_p = rho_0*g[2]*(x[2] - L[2]) 
        phi = bubbleRadius - math.sqrt((x[0] - bubbleCenter_x)**2 + (x[1] - bubbleCenter_y)**2+(x[2] - bubbleCenter_z)**2)
        if phi < 0.0:
            return pure_water_p
        if phi >= 0.0:
            return pure_water_p+sigma_01/bubbleRadius

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

#level set

class Bubble_phi:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        import math
        dx = ellipticalBubble_x*(X[0]-self.center[0]);
        dy = ellipticalBubble_y*(X[1]-self.center[1]);
        dz = ellipticalBubble_z*(X[2]-self.center[2])
        dBubble = self.radius - math.sqrt(dx**2 + dy**2 + dz**2)
        dwaterLevel = X[2] - waterLevel
        if math.fabs(dBubble) < math.fabs(dwaterLevel):
            return dBubble
        else:
            return dwaterLevel
    def uOfXT(self,X,t):
        return self.uOfX(X)

def getDBC_phi(x):
    pass

def getAFBC_phi(x):
     pass
 
#vof

def getDBC_vof(x):
    pass

def getAFBC_vof(x):
    pass

def Heaviside(phi):
    if phi > 0:
        return 1.0
    else:
        return 0.0

class Bubble_H:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.bubble_phi = Bubble_phi(center,radius)
    def uOfX(self,X):
        return Heaviside(self.bubble_phi.uOfX(X))
    def uOfXT(self,X,t):
        return self.uOfX(X)
