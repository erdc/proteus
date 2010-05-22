from math import sin,cos,pi,fabs,exp,sqrt
from pyadh.ctransportCoefficients import smoothedHeaviside
import wavetank2d_gen_poly
import pyadh
applyCorrection=True
applyRedistancing=True
nd = 2
wavetank_quad_order = 3
runCFL = 0.33

useStokes=False
nLevels = 1
epsFact = 1.5
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_curvature = 1.5
epsFact_redistance = 0.33
epsFact_massCorrection_heaviside = 1.5
epsFact_massCorrection_dirac = 1.5
epsFact_massCorrection_diffusion = 10.0

ns_lagSubgridError=True
ns_shockCapturingFactor=0.33
ls_shockCapturingFactor=0.33
rd_shockCapturingFactor=0.33
vof_shockCapturingFactor=0.33
#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1=1.500e-5
#surface tension
sigma_01 = 0.0#72.8e-3
#gravity
g=[0.0,-9.8]

polyfile = "wavetank"

tankHeight=15.0
beachSlope = 5.0/100.0
beachStartX = 1.0+10.0+25.0
boundaryTags = wavetank2d_gen_poly.genPoly(polyfile,
                                           water_sponge_length=2.0,
                                           piston_length=20.0,
                                           bottom_length=25.0,
                                           beach_slope=beachSlope,
                                           beach_length_x=100.0,
                                           land_length=beachStartX/3.0,
                                           land_sponge_length=beachStartX/3.0,
                                           back_height=10.0)

he = 0.02*tankHeight
triangleOptions="q30Dena%f" % (0.5*he**2)

waterLevel = 5.0#0.5*tankHeight
waveAmplitude = 3.0
wavePeriod = 10.0 #s
dt_init = 0.1*wavePeriod
generateWaves=True
outFlowEnd = False#True
#Navier-Stokes

T=2*wavePeriod

#inflow flux and  velocity for wave generator boundary
def waveGenerator_v(x,t):
    if generateWaves:
        return waveAmplitude*(2.0*pi/wavePeriod)*sin(2.0*pi*t/wavePeriod)#*cos((pi/2.0)*(x[0]/pistonEnd)) #smooth down to no flow 
    else:
        return 0.0

def waveGenerator_flux(x,t):
    if generateWaves:
        return -waveAmplitude*(2.0*pi/wavePeriod)*sin(2.0*pi*t/wavePeriod)#*cos((pi/2.0)*(x[0]/pistonEnd)) #smooth down to no flow 
    else:
        return 0.0

EPS = 1.0e-8

def getDBC_p_wavetank(x,flag):
    if flag==boundaryTags['top']:
        return lambda x,t: 0.0

def getDBC_u_wavetank(x,flag):
    if flag == boundaryTags['piston']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['walls']:
        return lambda x,t: 0.0

def getDBC_v_wavetank(x,flag):
    if flag == boundaryTags['piston']:
        return waveGenerator_v
    elif flag == boundaryTags['walls']:
        return lambda x,t: 0.0

def getDBC_vof_wavetank(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 1.0
    elif flag == boundaryTags['walls']:
        return lambda x,t: 0.0

def getAFBC_p_wavetank(x,flag):
    if flag == boundaryTags['piston']:
        return waveGenerator_flux
    elif flag == boundaryTags['walls']:
       return lambda x,t: 0.0

def getAFBC_vof_wavetank(x,flag):
    if flag == boundaryTags['top']:
        return None
    else:
        lambda x,t: 0.0

def getAFBC_u_wavetank(x,flag):
    pass
#    if flag == boundaryTags['walls']:
#        return lambda x,t: 0.0

def getAFBC_v_wavetank(x,flag):
    pass
#    if flag == boundaryTags['walls']:
#        return lambda x,t: 0.0

def getDFBC_u_wavetank(x,flag):
   if flag == boundaryTags['top']:
       return lambda x,t: 0.0
#    if flag in [boundaryTags['walls'],boundaryTags['top']]:
#        return lambda x,t: 0.0
    
def getDFBC_v_wavetank(x,flag):
   if flag == boundaryTags['top']:
       return lambda x,t: 0.0
#    if flag in [boundaryTags['walls'],boundaryTags['top']]:
#        return lambda x,t: 0.0

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= waterLevel:
            return -(tankHeight-x[1])*rho_1*g[1]
        else:
            return -((tankHeight-waterLevel)*rho_1 +
                     (waterLevel-x[1])*rho_0)*g[1]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

#level set
def getDBC_phi(x,flag):
    pass

class Flat_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        wby=waterLevel
        wbx=beachStartX +  waterLevel/beachSlope
        if x[0] > wbx:
            return sqrt((wbx-x[0])**2 + (wby-x[1])**2)
        else:
            return x[1] - waterLevel

class Flat_H:
    def __init__(self):
        self.phi=Flat_phi()
    def uOfXT(self,x,t):
        #return (1.0-smoothedHeaviside(he*epsFact_massCorrection_heaviside,self.phi.uOfXT(x,t)))
        return smoothedHeaviside(he*epsFact_massCorrection_heaviside,self.phi.uOfXT(x,t))

restrictFineSolutionToAllMeshes=False
parallelPartitioningType = pyadh.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 2
