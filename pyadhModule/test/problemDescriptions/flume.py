from math import exp
from pyadh.ctransportCoefficients import smoothedHeaviside

applyCorrection=True
applyRedistancing=True
nd = 2


#L = (1.0,1.0,1.0)
#L = (1.0e-1,1.0e-1,1.0)
#L = (1.0e-2,1.0e-2,1.0)
#L = (1.0e-3,1.0e-3,1.0)
#L = (2.0,2.0,1.0)
#polyfile="pome_haines_cube"
#L = (1.0e-2,3.0*1.0e-2)

polyfile = "flume"
nLevels=1
flume_quad_order = 3
runCFL = 1.0
nDTout = 100

ns_lagSubgridError = True
ns_shockCapturingFactor = 0.33
ls_shockCapturingFactor = 0.33
lag_ls_shockCapturing = True
vof_shockCapturingFactor = 0.33
rd_shockCapturingFactor = 0.33

usePETSc=True

useStokes=False
epsFact = 1.5
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_curvature = 1.5
epsFact_redistance = 0.33
epsFact_massCorrection_heaviside = 1.5
epsFact_massCorrection_dirac = 1.5
epsFact_massCorrection_diffusion = 10.0
#epsFact = 0.0
#water
#sigma_01 = 72.8e-3
sigma_01 = 0.0
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1=1.500e-5
#
#rho_1 = rho_0
#rho_1 = rho_0
#nu_1 = nu_0
#gravity
g=[0.0,-9.8]
#g=[0.0,0.0]


f = open(polyfile+'.poly','r')
lines = f.readlines()
f.close()
in_grate_right = float(lines[3].split()[1])
weirHeight = float(lines[5].split()[2])
weirStart = float(lines[5].split()[1])
weirEnd = float(lines[6].split()[1])
out_grate_left = float(lines[8].split()[1])
flumeEnd = float(lines[9].split()[1])
flumeTop = float(lines[10].split()[2])

#he = 0.065*flumeTop
he = 0.1*flumeTop
#he = 0.25*flumeTop
triangleOptions="q30Dena%f" % (0.5*he**2,)
bf = {'inflow_grate':1,
                 'outflow_grate':4,
                 'outflow':5,
                 'bottom':2,
                 'weir':3,
                 'top':6,
                 'inflow':7}
waterLevel = 0.5*weirHeight + 0.0*(flumeTop - weirHeight)
waterLevel_init = 1.5*weirHeight
inflowStop = 0.5*flumeTop

average_u = 0.1
#inflow = average_u*waterLevel/in_grate_right #
#rampDT=10.0
T = 20.0*flumeEnd/average_u
dt_init = 0.1*flumeEnd/average_u

if average_u > 0.0:
    rampDT=flumeEnd/average_u
def rampInflow(t):
    return 1.0 - exp(-2.4 * t/rampDT)

print in_grate_right,weirHeight,weirStart,weirEnd,flumeEnd,flumeTop,out_grate_left,triangleOptions


EPS = 1.0e-8

def onWeir(x):
    if (x[0] >= weirStart-EPS and
        x[0] <= weirEnd+EPS and
        x[1] <= weirHeight + EPS):
        return True
    else:
        return False

def onBottom(x):
    if(x[1] <= EPS):
        return True
    else:
        return False

def onTop(x):
    if (x[1] >= flumeTop - EPS):
        return True
    else:
        return False

def onSides(x):
    if (x[0] <= EPS or
        x[0] >= flumeEnd - EPS):
        return True
    else:
        return False

def onInflow(x):
    if onBottom(x):
        if (x[0] >= 0.0 + EPS
            and
            x[0] <= in_grate_right - EPS):
            return True
        else:
            return False
    else:
        return False
    
def onOutflow(x):
    if onBottom(x):
        if (x[0] >= out_grate_left + EPS
            and
            x[0] <= flumeEnd - EPS):
            return True
        else:
            return False
    else:
        return False

def getDBC_p_weir(x,flag):
    if flag == bf['top']:
        return lambda x,t: 0.0
    elif flag == bf['outflow']:
        return lambda x,t: -(flumeTop-x[1])*rho_1*g[1]
#    if onTop(x): #top open
#        return lambda x,t: 0.0

walls = [bf['bottom'],bf['weir'],bf['inflow_grate'],bf['outflow_grate']]
def getDBC_u_weir(x,flag):
    if flag in walls:
        return lambda x,t: 0.0
#    if (onBottom(x) or onWeir(x) or onSides(x)):
#        return lambda x,t: 0.0

def getDBC_v_weir(x,flag):
    if flag in walls:
        return lambda x,t: 0.0
#    if onSides(x):
#        return lambda x,t: 0.0
#    if (onBottom(x) or onWeir(x)):
#        if onInflow(x):
#            return lambda x,t: inflow*rampInflow(t)
#        elif onOutflow(x):
#            return lambda x,t: -inflow*rampInflow(t)
#            #return lambda x,t: inflow*rampInflow(t) #cek debug
#        else:
#            return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x,flag):
    if flag == bf['inflow']:
        return lambda x,t: -average_u#inflow
    elif flag in walls:
        return lambda x,t: 0.0
#   if onTop(x):
#       pass
#   elif onBottom(x):
#       if onInflow(x):
#           return lambda x,t: -inflow*rampInflow(t)
#       elif onOutflow(x):
#           return lambda x,t: inflow*rampInflow(t)
#           #return lambda x,t: -inflow*rampInflow(t) #cek debug
#       else:
#           return lambda x,t: 0.0
#   else:
#        return lambda x,t: 0.0

def getAFBC_u_weir(x,flag):
    pass

def getAFBC_v_weir(x,flag):
    pass

def getDFBC_u_weir(x,flag):
    if flag == bf['top']:
        return lambda x,t: 0.0
#    if onTop(x):
#        return lambda x,t: 0.0

def getDFBC_v_weir(x,flag):
    if flag == bf['top']:
        return lambda x,t: 0.0
#    if onTop(x):
#        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir,
                                    1:getAFBC_u_weir,
                                    2:getAFBC_v_weir}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_weir,2:getDFBC_u_weir},
                                   2:{1:getDFBC_v_weir,2:getDFBC_v_weir}}

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= waterLevel:
            return -(flumeTop-x[1])*rho_1*g[1]
        else:
            return -((flumeTop-waterLevel)*rho_1 +
                     (waterLevel-x[1])*rho_0)*g[1]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
    
class NoWeir_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0#average_u

class NoWeir_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Hydrostatic_p(),
                     1:NoWeir_u(),
                     2:NoWeir_v()}

#level set
def getDBC_phi(x,flag):
    pass

class Flat_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return x[1] - waterLevel_init

#vof
def Heaviside(phi):
    if phi > 0:
        return 1.0
    else:
        return 0.0

class Flat_H:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_massCorrection_heaviside*he,(x[1] - waterLevel_init))

def getDBC_vof(x,flag):
    if flag == bf['inflow']:
        if x[1] < waterLevel:
            return lambda x,t: 0.0
        else:
            return lambda x,t: 1.0
    elif flag == bf['outflow']:
        return lambda x,t: 1.0
#    if onTop(x):
#        return lambda x,t: 1.0
#    else:
#        return None

def getAFBC_vof(x,flag):
    if flag in walls:
        return lambda x,t: 0.0
#    if onTop(x):
##        return None
#    else:
#        return lambda x,t: 0.0
