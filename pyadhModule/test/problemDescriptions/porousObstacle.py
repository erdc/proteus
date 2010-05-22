from math import sqrt,pi,sin,cos
applyCorrection=True
applyRedistancing=True
nd = 2

#flume 1
U0     = 0.1  #m/s
height = 1.2  #m
length = 51.0 #m
length2height = length/height#10.0

in_grate_right = 1.0
out_grate_left = 50.0

waterLevel = 0.8

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

Re = height*U0/nu_0
Fr = U0/sqrt(height*gms2)

meanGrainsize = 0.05 #m
obstaclePorosity = 0.5#0.45


polyfile = 'cobras_wave_flume_1'

porosityTypes = {0 : 1.0, #fluid
                 1 : 0.5, #tetrapod
                 2 : 0.53,#filter_layer
                 3 : 0.49,#rubble_mound
                 4 : 1.0e-12}#caisson not meshed
meanGrainTypes= {0 : 0.05,#fluid shouldn't matter
                 1 : 0.05, #tetrapod
                 2 : 0.05,#filter_layer
                 3 : 0.05,#rubble_mound
                 4 : 0.05}#caisson not meshed




average_u = U0
inflow    = average_u*waterLevel/in_grate_right


nn=3
nLevels=2 #1 very coarse
porousObstacle_quad_order = 4
runCFL = 10.0
DT = 1.0e-2
T=9.0*DT#9999*DT
nDTout = int(round(T/DT))

epsFact = 1.5
epsFact_density_viscosity = 1.5
epsFact_curvature = 3.0
epsFact_redistance = 0.5#smoothing in right hand side of eikonal equation (switches from 1 to -1)
epsFact_massCorrection_heaviside = 1.0e-2
epsFact_massCorrection_dirac = 1.0e-2
epsFact_massCorrection_diffusion = 1.0e2 #not a typo


EPS = 1.0e-8
#from cobras_wave_flume_1
segmentLabels = {'default':0,
                 'left': 1,
                 'bottom' : 2,
                 'right'  : 3,
                 'top'    : 4,
                 'bottom_slope': 5,
                 'bottom_shelf': 6,
                 'tetrapod_exterior':7,
                 'caisson_exterior' :8,
                 'tetrapod_interior':9,
                 'caisson_interior' :9,
                 'rubble_mound': 10,
                 'filter_layer_exterior': 11,
                 'bottom_inflow' : 12,
                 'bottom_outflow' : 13}

#collect some tags
bottomBoundaries = [segmentLabels['bottom'],segmentLabels['bottom_slope'],segmentLabels['bottom_shelf'],
                    segmentLabels['bottom_inflow'],segmentLabels['bottom_outflow']]
                    

#use 
#weir (constant flow bcs)
#or wave conditions
useWaveConditions = True
generateWaves = True
waveAmplitude = 0.33*waterLevel
wavePeriod = 1.0/5.0
pistonEnd  = in_grate_right
def getDBC_p_weir(x,tag):
    if tag == segmentLables['top']:
        return lambda x,t: 0.0

def getDBC_u_weir(x,tag):
    if (tag in bottomBoundaries or
        tag in [segmentLabels['left'],segmentLabels['right'],
                segmentLabels['caisson_exterior']]):
        return lambda x,t: 0.0
    #if x[1] <= EPS:
    #    if x[0] <= in_grate_right:
    #        return lambda x,t: 0.0
    #    elif x[0] >= out_grate_left:
    #        return lambda x,t: 0.0
    #make left and right hand boundaries slip?
    #if x[0] <= EPS:
    #    return lambda x,t: 0.0
    #if x[0] >= length-EPS:
    #    return lambda x,t: 0.0

def getDBC_v_weir(x,tag):
    if tag == segmentLabels['bottom_inflow']:
        return lambda x,t: inflow
    elif tag == segmentLabels['bottom_outflow']:
        return lambda x,t: -inflow
    elif tag in [segmentLables['bottom'],segmentLabels['bottom_slope'],segmentLabels['bottom_shelf']]:
        return lambda x,t: 0.0
        
    #if x[1] <= EPS:
    #    if x[0] <= in_grate_right:
    #        return lambda x,t: inflow
    #    elif x[0] >= out_grate_left:
    #        return lambda x,t: -inflow
    #    else: #
    #        return lambda x,t: 0.0
    
#         else:
#             return lambda x,t: 0.0 #no slip
#     elif not (x[1] >= height - EPS): #no slip everywhere else
#         return lambda x,t: 0.0

def getAFBC_p_weir(x,tag):
    if tag == segmentLabels['bottom_inflow']:
        return lambda x,t: -inflow
    elif tag == segmentLabels['bottom_outflow']:
        return lambda x,t: inflow
    elif tag in [segmentLables['bottom'],segmentLabels['bottom_slope'],segmentLabels['bottom_shelf']]:
        return lambda x,t: 0.0
    elif tag == segmentLabels['top']:
        pass
    else: #check this
        return lambda x,t: 0.0
    #if x[1] <= EPS:
    #    if x[0] <= in_grate_right:
    #        return lambda x,t: -inflow
    #    elif x[0] >= out_grate_left:
    #        return lambda x,t: inflow
    #    else:
    #        return lambda x,t: 0.0
    #elif x[1] >= height - EPS:
    #    pass
    #else: #slip elsewhere
    #    return lambda x,t: 0.0


#generate some waves

#inflow flux and  velocity for wave generator boundary
def waveGenerator_v(x,t):
    if generateWaves:
        return waveAmplitude*2.0*pi*wavePeriod*sin(2.0*pi*t*wavePeriod)*cos((pi/2.0)*(x[0]/pistonEnd)) #smooth down to no flow 
    else:
        return 0.0

def waveGenerator_flux(x,t):
    if generateWaves:
        return -waveAmplitude*2.0*pi*wavePeriod*sin(2.0*pi*t*wavePeriod)*cos((pi/2.0)*(x[0]/pistonEnd)) #smooth down to no flow 
    else:
        return 0.0

EPS = 1.0e-8


def onBottom(x,tag):
    if tag in bottomBoundaries:
        return True
    else:
        return False

def onPiston(x,tag):
    if tag == segmentLabels['bottom_inflow']:
        return True
    else:
        return False

def onSides(x,tag):
    if tag in [segmentLabels['left'],segmentLabels['right']]: 
        return True
    else:
        return False
    
def onTop(x,tag):
    if tag == segmentLabels['top']:
        return True
    else:
        return False

def getDBC_p_wavetank(x,tag):
    if onTop(x,tag):
        return lambda x,t: 0.0

def getDBC_u_wavetank(x,tag):
    if onTop(x,tag):
        return lambda x,t: 0.0
#     if onBottom(x) or onSides(x): #no slip
#         return lambda x,t: 0.0

def getDBC_v_wavetank(x,tag):
    if onPiston(x,tag):
        return waveGenerator_v
#     elif onBottom(x):# or onSides(x): #no slip
#         return lambda x,t: 0.0

def getAFBC_p_wavetank(x,tag):
    if onPiston(x,tag):
        return waveGenerator_flux
    elif (onBottom(x,tag) or onSides(x,tag)):
        return lambda x,t: 0.0
    else:
        pass

def getAFBC_u_wavetank(x,tag):
    pass

def getAFBC_v_wavetank(x,tag):
    pass


if useWaveConditions:
    getDBC_p = getDBC_p_wavetank
    getDBC_u = getDBC_u_wavetank
    getDBC_v = getDBC_v_wavetank
    getAFBC_p= getAFBC_p_wavetank
else:
    getDBC_p = getDBC_p_weir
    getDBC_u = getDBC_u_weir
    getDBC_v = getDBC_v_weir
    getAFBC_p= getAFBC_p_weir
 

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p}

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
def getDBC_phi(x,tag):
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

def getDBC_vof(x,tag):
    pass

def getAFBC_vof(x,tag):
    pass
