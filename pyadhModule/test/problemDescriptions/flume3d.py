from math import exp
applyCorrection=True
applyRedistancing=True
nd = 3


polyfile = "flume3d"
nLevels=1
flume_quad_order = 3
runCFL = 1
nDTout = 100

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
sigma_01 = 0.0
#sigma_01 = 0.0
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
g=[0.0,0.0,-9.8]


f = open(polyfile+'.poly','r')
lines = f.readlines()
f.close()
in_grate_right = float(lines[3].split()[1])
weirHeight = float(lines[5].split()[3])
weirStart = float(lines[5].split()[1])
weirEnd = float(lines[6].split()[1])
out_grate_left = float(lines[8].split()[1])
flumeEnd = float(lines[9].split()[1])
flumeTop = float(lines[10].split()[3])
flumeBack = float(lines[11].split()[2])

he = 0.1*flumeTop
triangleOptions="q30Dena%f" % ((1.0/6.0)*(0.25*flumeTop)**3,)

waterLevel = weirHeight + 0.0*(flumeTop - weirHeight)
inflowStop = 0.5*flumeTop

average_u = 1.0
inflow = average_u*waterLevel/in_grate_right #
rampDT=10.0
if average_u > 0.0:
    rampDT=flumeEnd/average_u
def rampInflow(t):
    return 1.0 - exp(-2.4 * t/rampDT)

print in_grate_right,weirHeight,weirStart,weirEnd,flumeEnd,flumeTop,out_grate_left,triangleOptions


EPS = 1.0e-8

def onWeir(x):
    if (x[0] >= weirStart-EPS and
        x[0] <= weirEnd+EPS and
        x[2] <= weirHeight + EPS):
        return True
    else:
        return False

def onBottom(x):
    if(x[2] <= EPS):
        return True
    else:
        return False

def onTop(x):
    if (x[2] >= flumeTop - EPS):
        return True
    else:
        return False

def onSides(x):
    if (x[0] <= EPS or
        x[0] >= flumeEnd - EPS):
        return True
    if (x[1] <= EPS or
        x[1] >= flumeBack - EPS):
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

def getDBC_p_weir(x):
    if onTop(x): #top open
        return lambda x,t: 0.0

def getDBC_u_weir(x):
    if (onBottom(x) or onWeir(x) or onSides(x) or onTop(x)):
        return lambda x,t: 0.0
def getDBC_v_weir(x):
    if (onBottom(x) or onWeir(x) or onSides(x) or onTop(x)):
        return lambda x,t: 0.0

def getDBC_w_weir(x):
    if onSides(x):
        return lambda x,t: 0.0
    if (onBottom(x) or onWeir(x)):
        if onInflow(x):
            return lambda x,t: inflow*rampInflow(t)
        elif onOutflow(x):
            return lambda x,t: -inflow*rampInflow(t)
            #return lambda x,t: inflow*rampInflow(t) #cek debug
        else:
            return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir,
                       3:getDBC_w_weir}

def getAFBC_p_weir(x):
   if onTop(x):
       pass
   elif onBottom(x):
       if onInflow(x):
           return lambda x,t: -inflow*rampInflow(t)
       elif onOutflow(x):
           return lambda x,t: inflow*rampInflow(t)
       else:
           return lambda x,t: 0.0
   else:
        return lambda x,t: 0.0

def getAFBC_u_weir(x):
    pass

def getAFBC_v_weir(x):
    pass

def getAFBC_w_weir(x):
    pass

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir,
                                    1:getAFBC_u_weir,
                                    2:getAFBC_v_weir,
                                    3:getAFBC_w_weir}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{},
                                   2:{},
                                   3:{}}

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= waterLevel:
            return -(flumeTop-x[2])*rho_1*g[2]
        else:
            return -((flumeTop-waterLevel)*rho_1 +
                     (waterLevel-x[2])*rho_0)*g[2]

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

class NoWeir_w:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Hydrostatic_p(),
                     1:NoWeir_u(),
                     2:NoWeir_v(),
                     3:NoWeir_w()}

#level set
def getDBC_phi(x):
    pass

class Flat_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return x[2] - waterLevel

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
        return Heaviside(x[2] - waterLevel)

def getDBC_vof(x):
    pass

def getAFBC_vof(x):
    pass
