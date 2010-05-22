from math import exp
applyCorrection=True
applyRedistancing=True
nd = 2

U0     = 1.0  #m/s
height = 1.0 #m
length2height = 10.#10.0
length = length2height*height

L = (length2height,height/height,1.0)

nn=41
nLevels=1
lwi_dike_quad_order = 3
runCFL = 10000
DT = 0.05#0.1,0.05
T=50*DT#50
nDTout = int(round(T/DT))

useStokes=False
epsFact = 3.0
epsFact_density = 3.0
epsFact_viscosity = 3.0
epsFact_curvature = 0.001
epsFact_redistance = 0.5#smoothing in right hand side of eikonal equation (switches from 1 to -1)
epsFact_massCorrection_heaviside = 0.1
epsFact_massCorrection_dirac = 0.1
epsFact_massCorrection_diffusion = 10.0
#epsFact = 0.0
#water
sigma_01 = 0.0#72.8e-3
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
g=[0.0,-9.8]
#g=[0.0,0.0]

in_grate_right = 0.1*length2height
out_grate_left = 0.9*length2height

#print "NS porous Obstacle U0= %s Re= %s Fr= %s g= %s" % (U0,Re,Fr,g) 
meanGrainsize = 0.05 #0.05 #m
obstaclePorosity = 0.45#0.2#0.45#0.8#
useSquare = False

obstacleXM = 0.45*length2height
obstacleXP = 0.55*length2height
obstacleYM = 0.0
obstacleYP = 0.4


#trying LWI dike
if not useSquare:
    obstacleXM = 0.55*length2height
    obstacleXP = 0.6*length2height
slopeM = 1.0/6.0 #front slope
slopeP = 1.0/3.0 #back slope
obstacleDXM = (obstacleYP-obstacleYM)/slopeM
obstacleDXP = (obstacleYP-obstacleYM)/slopeP



waterLevel = obstacleYP*1.2
#waterLevel = obstacleYP*0.4
inflowStop = 0.5*L[1]

average_u = U0
inflow = average_u*waterLevel/in_grate_right #
rampDT=1.0
if average_u > 0.0:
    rampDT=L[0]/average_u
def rampInflow(t):
    return 1.0 - exp(-2.4 * t/rampDT)

#print in_grate_right,weirHeight,weirStart,weirEnd,flumeEnd,flumeTop,out_grate_left,triangleOptions


EPS = 1.0e-6

def onSquareObstacle(x):
    if obstacleXM <= x[0] and x[0] <= obstacleXP:
        if obstacleYM <= x[1] and x[1] <= obstacleYP:
            return True
    return False

def onLWIdike(x):
    if obstacleXM <= x[0] and x[0] <= obstacleXP:
        if obstacleYM <= x[1] and x[1] <= obstacleYP:
            return True
    if obstacleXM-obstacleDXM <= x[0] and x[0] <= obstacleXM:
        if obstacleYM <= x[1] and x[1] <= obstacleYM+slopeM*(x[0]-(obstacleXM-obstacleDXM)):
            return True
    if obstacleXP <= x[0] and x[0] <= obstacleXP+obstacleDXP:
        if obstacleYM <= x[1] and x[1] <= obstacleYP+slopeP*(obstacleXP-x[0]):
            return True
    return False

def setObstaclePorosity(x,porosity):
    porosity.flat[:]  = 1.0 #fluid domain
    #import pdb
    #pdb.set_trace()
    nPoints = len(x.flat)/3 #points always 3d
    for k in range(nPoints):
        if useSquare and onSquareObstacle(x.flat[3*k:3*(k+1)]):
            porosity.flat[k]=obstaclePorosity
        elif onLWIdike(x.flat[3*k:3*(k+1)]):
            porosity.flat[k]=obstaclePorosity
        #
    #
    #mwf hack
    #porosity.flat[:] = obstaclePorosity
#


def onBottom(x):
    if(x[1] <= EPS):
        return True
    else:
        return False

def onTop(x):
    if (x[1] >= L[1] - EPS):
        return True
    else:
        return False

def onSides(x):
    if (x[0] <= EPS or
        x[0] >= L[0] - EPS):
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
            x[0] <= L[0] - EPS):
            return True
        else:
            return False
    else:
        return False

def getDBC_p_lwi_dike(x):
    if onTop(x): #top open
        return lambda x,t: 0.0

def getDBC_u_lwi_dike(x):
    if (onBottom(x) or onSides(x) or onTop(x)):
        return lambda x,t: 0.0

def getDBC_v_lwi_dike(x):
    if onSides(x):
        return lambda x,t: 0.0
    if (onBottom(x)):
        if onInflow(x):
            return lambda x,t: inflow*rampInflow(t)
        elif onOutflow(x):
            return lambda x,t: -inflow*rampInflow(t)
            #return lambda x,t: inflow*rampInflow(t) #cek debug
        else:
            return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_lwi_dike,
                       1:getDBC_u_lwi_dike,
                       2:getDBC_v_lwi_dike}

def getAFBC_p_lwi_dike(x):
   if onTop(x):
       pass
   elif onBottom(x):
       if onInflow(x):
           return lambda x,t: -inflow*rampInflow(t)
       elif onOutflow(x):
           return lambda x,t: inflow*rampInflow(t)
           #return lambda x,t: -inflow*rampInflow(t) #cek debug
       else:
           return lambda x,t: 0.0
   else:
        return lambda x,t: 0.0

def getAFBC_u_lwi_dike(x):
    pass

def getAFBC_v_lwi_dike(x):
    pass

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_lwi_dike,1:getAFBC_u_lwi_dike,1:getAFBC_v_lwi_dike}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= waterLevel:
            return -(L[1]-x[1])*rho_1*g[1]
        else:
            return -((L[1]-waterLevel)*rho_1 +
                     (waterLevel-x[1])*rho_0)*g[1]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
    
class NoDike_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0#average_u

class NoDike_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Hydrostatic_p(),
                     1:NoDike_u(),
                     2:NoDike_v()}

#level set
def getDBC_phi(x):
    pass

class Flat_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return x[1] - waterLevel

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
        return Heaviside(x[1] - waterLevel)

def getDBC_vof(x):
    pass

def getAFBC_vof(x):
    pass
