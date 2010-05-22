from math import sqrt
applyCorrection=True
applyRedistancing=True
nd = 2

height = 4.0  #m 4.0
length = 20.0 #m
L = (length,height,1.0)
obstacleXM = 0.45*length
obstacleXP = 0.55*length
obstacleYM = 0.0
obstacleYP = min(1.0,0.5*height)
in_grate_right = 1.0
out_grate_left = 19.0

waterLevel = obstacleYP + 0.05*(height-obstacleYP)#0.5*(height-obstacleYP)
inflowStop = 0.5*height

#epsFact = 0.0
#water
sigma_01 = 72.8e-3
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
gms2 = 9.8
g=[0.0,-gms2]
#g=[0.0,0.0]
meanGrainsize = 0.05 #m
obstaclePorosity = 0.15#0.45

def setObstaclePorosity(x,porosity,meanGrain=None):
    porosity.flat[:]  = 1.0 #fluid domain
    if meanGrain != None:
        meanGrain.flat[:] = meanGrainsize
    #import pdb
    #pdb.set_trace()
    nPoints = len(x.flat)/3 #points always 3d
    for k in range(nPoints):
        if obstacleXM <= x.flat[k*3+0] and x.flat[k*3+0] <= obstacleXP:
            if obstacleYM <= x.flat[k*3+1] and x.flat[k*3+1] <= obstacleYP:
                porosity.flat[k]=obstaclePorosity
            #
        #
    #
    #mwf hack
    #porosity.flat[:] = obstaclePorosity
#

U0     = 0.1 #m/s

Re = height*U0/nu_0
Fr = U0/sqrt(height*gms2)


average_u = U0
inflow    = average_u*waterLevel/in_grate_right


nn=41
nLevels=1
porousObstacle_quad_order = 4
runCFL = 10.0
DT = 1.0e-2
T=9999*DT
nDTout = int(round(T/DT))

epsFact = 1.5
epsFact_density_viscosity = 1.5
epsFact_curvature = 3.0
epsFact_redistance = 0.5#smoothing in right hand side of eikonal equation (switches from 1 to -1)
epsFact_massCorrection_heaviside = 1.0e-2
epsFact_massCorrection_dirac = 1.0e-2
epsFact_massCorrection_diffusion = 1.0e2 #not a typo

#f = open(polyfile+'.poly','r')
#lines = f.readlines()
#f.close()
#in_grate_right = float(lines[3].split()[1])
#weirHeight = float(lines[5].split()[2])
#weirStart = float(lines[5].split()[1])
#weirEnd = float(lines[6].split()[1])
#out_grate_left = float(lines[8].split()[1])
#flumeEnd = float(lines[9].split()[1])
#flumeTop = float(lines[10].split()[2])



EPS = 1.0e-8

def getDBC_p_weir(x):
    if x[1] >= height - EPS:
        return lambda x,t: 0.0

def getDBC_u_weir(x):
    if x[1] <= EPS:
        if x[0] <= in_grate_right:
            return lambda x,t: 0.0
        elif x[0] >= out_grate_left:
            return lambda x,t: 0.0
    #make left and right hand boundaries slip?
    if x[0] <= EPS:
        return lambda x,t: 0.0
    if x[0] >= length-EPS:
        return lambda x,t: 0.0
    #elif x[1] >= height - EPS:
    #    return lambda x,t: 0.0
#         else:
#             return lambda x,t: 0.0 #no slip
#     elif not (x[1] >= height - EPS): #no slip everywhere else
#         return lambda x,t: 0.0

def getDBC_v_weir(x):
    if x[1] <= EPS:
        if x[0] <= in_grate_right:
            return lambda x,t: inflow
        elif x[0] >= out_grate_left:
            return lambda x,t: -inflow
        else: #enforce slip on bottom
            return lambda x,t: 0.0
    
#         else:
#             return lambda x,t: 0.0 #no slip
#     elif not (x[1] >= height - EPS): #no slip everywhere else
#         return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x):
    if x[1] <= EPS:
        if x[0] <= in_grate_right:
            return lambda x,t: -inflow
        elif x[0] >= out_grate_left:
            return lambda x,t: inflow
        else:
            return lambda x,t: 0.0
    elif x[1] >= height - EPS:
        pass
    else: #slip elsewhere
        return lambda x,t: 0.0
#     if x[1] <= EPS:
#         if x[0] <= inflowStop:
#             return lambda x,t: -inflow
#     elif x[0] >= flumeEnd - EPS:
#         pass
#     elif x[1] >= height - EPS:
#         pass
#     else:
#         return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir}

diffusiveFluxBoundaryConditions = {0:{}}

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= waterLevel:
            return -(height-x[1])*rho_1*g[1]
        else:
            return -((height-waterLevel)*rho_1 +
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
        return average_u

class NoWeir_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Hydrostatic_p(),
                     1:NoWeir_u(),
                     2:NoWeir_v()}

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
