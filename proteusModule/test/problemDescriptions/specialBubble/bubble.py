applyCorrection=True
applyRedistancing=True
lag_ns_subgridError=True
rd_shockCapturingFactor=0.9
nd = 2
L = (1.0,1.0,1.0)

#nnx=81
#nny=161
nnx=21
nny=21
nLevels=1
spaceOrder=1
bubble_quad_order = 3
#runCFL = 100000.0
#T=0.5
T=100.0
#T=6.0e-3
#T=0.001
#DT=1.0e-2
nDTout = 1001
nDTout = None
useBackwardEuler=True
runCFL=0.33
usePETSc=True

useStokes=False
epsFact_density = 3.0
epsFact_viscosity =3.0
rdtimeIntegration = 'newton'
epsFact_redistance = 0.33
rdtimeIntegration = 'osher'
epsFact_curvature = 0.33
epsFact_consrv_heaviside = 3.0
epsFact_consrv_dirac = 3.0
epsFact_consrv_diffusion=3.0

#water
sigma_01 = 72.8e-3
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1=1.500e-5

#gravity
g=[0.0,0.0]

waterLevel =1.75*L[1]
import math
bubbleRadius = 0.25*L[0]
bubbleCenter_x = 0.5*L[0]
bubbleCenter_y = 2.0*bubbleRadius

def getDBC_p_bubble(x,flag):
    import math
    #top open
    if x[1] == L[1]:
        return lambda x,t: 0.0

EPS = 1.0e-8
def getDBC_u_bubble(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] == 0.0):
        return lambda x,t: 0.0
#    #slip
#    pass

def getDBC_v_bubble(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] == 0.0):
        return lambda x,t: 0.0
#    #slip
#    pass

dirichletConditions = {0:getDBC_p_bubble,
                       1:getDBC_u_bubble,
                       2:getDBC_v_bubble}

def getAFBC_p_bubble(x,flag):
    #bottom and sides no flow
    if (x[0] in [0.0,L[0]] or
        x[1] == 0.0):
        return lambda x,t: 0.0

def getAFBC_u_bubble(x,flag):
    #bottom and sides no flow
    if (x[0] in [0.0,L[0]] or
        x[1] == 0.0):
        return lambda x,t: 0.0

def getAFBC_v_bubble(x,flag):
    #bottom and sides no flow
    if (x[0] in [0.0,L[0]] or
        x[1] == 0.0):
        return lambda x,t: 0.0

def getDFBC_u_bubble(x,flag):
    return lambda x,t: 0.0

def getDFBC_v_bubble(x,flag):
    return lambda x,t: 0.0
    
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_bubble,
                                    1:getAFBC_u_bubble,
                                    2:getAFBC_v_bubble}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_bubble},
                                   2:{2:getDFBC_v_bubble}}

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        import math
        return 0.0
        phi = bubbleRadius - math.sqrt((x[0] - bubbleCenter_x)**2 + (x[1] - bubbleCenter_y)**2)
        if phi < 0.0:
            return 0.0
        if phi >= 0.0:
            return sigma_01/bubbleRadius

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Hydrostatic_p(),
                     1:AtRest(),
                     2:AtRest()}

