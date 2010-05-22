from math import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from pyadh import Domain
from pyadh import TransportCoefficients
from pyadh.TransportCoefficients import TwophaseNavierStokes_ST_LS_SO,VOFCoefficients,TC_base
import beach_erosion_board_2dDomain
#numerics options
applyCorrection  =True#True#False
applyRedistancing=True#right now the solution is more symmetric without redistancing
rdtimeIntegration='tte'
sloshbox_quad_order = 4 #need 4 for quadratics

useBackwardEuler=True#False
useBackwardEuler_ls=True#False

lag_ns_subgridError=True
lag_ls_shockCapturing=True
ls_shockCapturingFactor=0.33#0.1
rd_shockCapturingFactor=0.33#0.1
ns_shockCapturingFactor=0.33

rd_freezeLS = True
usePETSc=True#False#True
useStokes=False
VOF = False
nonconsrv = True
if  VOF:
    nonconsrv = True
spaceOrder=1

#spatial domain
nd = 2


height=1.0 #should be 0.2125
length= 9.98/2.0#9.98
L = (length,height) #[m]
domain = Domain.RectangularDomain(L=L,
                                  name="LinLuiEx1")
#grid resolution 
nnx = 41
nny = 41
nLevels=1
dy = domain.L[1]/float(nny-1.)
dx = domain.L[0]/float(nnx-1.)

domain.writeAsymptote("LinLiuDomainEx1")
domain.writePoly("LinLiuDomainEx1")



#wave specification
runWaveProblem=True

#source zone for generating waves

#background water level
waterLevelBase = 0.5*height#should be 0.2 #m
waveHeight=0.2*waterLevelBase#0.05*waterLevelBase #[m]
wavePeriod= 1.   #[s]
waveCelerity  = 1.21#[m/s]
waveNumber= 1.04/waterLevelBase #[1/m]
waveLength= 2.0*pi/waveNumber
waveFrequency = 2.0*pi/wavePeriod
source_height=max(5*dy,0.25*waterLevelBase)
source_ym    =domain.x[1]+0.3*waterLevelBase
Omega_s=[[L[0]*0.5-dx,L[0]*0.5+dx],
         [source_ym-0.5*source_height,source_ym+0.5*source_height],
         [0.0,1.0]]
   
sourceVolume = (Omega_s[0][1]-Omega_s[0][0])*(Omega_s[1][1]-Omega_s[1][0])
waveFlag= 0 #0 -- monochromatic
            #1 -- second order stokes
            #2 -- solitary wave
            #-1 -- no internal wave
setWavesAtInflow = False#True
outFlowEnd = True
if setWavesAtInflow:
    waveFlag = -1
height = domain.L[1]
length = domain.L[0]

topOpen = True#False


#boundary conditions
#slip on bottom, 
#outflow on sides
#top open with grate --> p = 0, u = 0, v is free
 

#
#numerical fudge factors
epsFact = 1.5
epsFact_density = 1.5#1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.33#1.5
epsFact_curvature=0.001
epsFact_source = 0.5
epsFact_vof=1.5
#get names right
if runWaveProblem:
    epsFact_redistance = 0.33
    epsFact_consrv_heaviside = 1.5
    epsFact_consrv_dirac = 1.5
    epsFact_consrv_diffusion = 1.5
else:
    epsFact_redistance = 1.5
    epsFact_consrv_heaviside=1.5
    epsFact_consrv_dirac=1.5
    epsFact_consrv_diffusion=1.5
    epsFact_vof=0.0


#physical parameters
#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5
#gravity
g=[0.0,-9.8]
#mwf play with density and viscosity ratios
#rho_1 = rho_0/20.
#nu_1  = nu_0*0.5
turbulenceClosureFlag = 0#Smagorinsky
smagorinskyConstant = 0.1

#surface tension
sigma_01 = 0.0#72.8e-3

T = wavePeriod*10.0#*10.0#10.#*10.0#2.5#3.0#10.0
import math
nDTout = math.ceil(T/wavePeriod)*10+1#51#101#None#2#None, can't be 1
runCFL = 0.33

dt_init = min(0.01*wavePeriod,T)


#set up absorbing boundary at left side using volume averaged ns formulation
spongeLayerWidth = max(2.0*dx,0.1*length)
useSpongeLayer = True#False
spongeGrainSize= 0.01#1.0
spongePorosity = 0.5
killNonlinearDragInSpongeLayer = True
def setSpongeLayer(x,porosity,meanGrain=None):
    porosity.fill(1.0) #fluid domain
    if meanGrain != None:
        meanGrain.fill(spongeGrainSize)
    #import pdb
    #pdb.set_trace()
    nPoints = len(x.flat)/3 #points always 3d
    for k in range(nPoints):
        if x.flat[k*3+0] <= spongeLayerWidth:
            porosity.flat[k]=spongePorosity
            #
        #
    #
    #mwf hack
    #porosity.flat[:] = obstaclePorosity
#
if useSpongeLayer == True:
    spongeLayerFunc = setSpongeLayer
else:
    spongeLayerFunc = None




#try to vary free surface height at boundary
epsForVOF=epsFact
def freeSurfaceVOF_bc(x,t):
    if x[0] <= 0.0:
        #mwf hack to debug time integration problems
        teval = t
        #waterLevel = waterLevelBase + waveHeight*sin(waveFrequency*teval)
        eta = 0.5*waveHeight*math.cos(waveNumber*(x[0]+0.25*waveLength)-waveFrequency*t)
        waterLevel = waterLevelBase + eta
        vof_eps = smoothedHeaviside(epsForVOF*dy,x[1] - waterLevel)
        #mwf debug
        #print "freeSurfaceVOF x[1]= %s t=%s waterLevel=%s vof=%s " % (x[1],t,waterLevel,vof_eps)
        return vof_eps
def waveVelocity_u(x,t):
    waterLevel =  waterLevelBase + waveHeight*sin(waveFrequency*t)
    vof_eps = smoothedHeaviside(epsForVOF*dy,x[1] - waterLevel)
    term1 =  0.5*waveHeight*abs(g[1])*wavePeriod*math.cosh(waveNumber*waterLevel)/(waveLength*math.cosh(waveNumber*waterLevelBase))
    term2 =  math.cos(waveNumber*(x[0]+0.25*waveLength)-waveFrequency*t)
    return (1.0-vof_eps)*term1*term2
def waveFlux_u(x,t):
    waterLevel =  waterLevelBase + waveHeight*sin(waveFrequency*t)
    vof_eps = smoothedHeaviside(epsForVOF*dy,x[1] - waterLevel)
    term1 =  0.5*waveHeight*abs(g[1])*wavePeriod*math.cosh(waveNumber*waterLevel)/(waveLength*math.cosh(waveNumber*waterLevelBase))
    term2 =  math.cos(waveNumber*(x[0]+0.25*waveLength)-waveFrequency*t)
    return -(1.0-vof_eps)*term1*term2
def waveVelocity_v(x,t):
    waterLevel =  waterLevelBase + waveHeight*sin(waveFrequency*t)
    vof_eps = smoothedHeaviside(epsForVOF*dy,x[1] - waterLevel)
    term1 =  0.5*waveHeight*abs(g[1])*wavePeriod*math.sinh(waveNumber*waterLevel)/(waveLength*math.cosh(waveNumber*waterLevelBase))
    term2 =  math.sin(waveNumber*(x[0]+0.25*waveLength)-waveFrequency*t)
    return (1.0-vof_eps)*term1*term2

def getDBC_p_sloshbox(x,flag):
    if x[1] == L[1]:
        return lambda x,t: 0.0

def getDBC_u_sloshbox(x,flag):
    if not outFlowEnd and x[0] >= L[0] - 1.0e-8:
        return lambda x,t: 0.0
    if setWavesAtInflow and x[0] <= 1.0e-8:
        return waveVelocity_u
    return None

def getDBC_v_sloshbox(x,flag):
    if x[1] <= domain.x[1]:
        return lambda x,t:0.0
    if setWavesAtInflow and x[0] <= 1.0e-8:
        return waveVelocity_v
def getDBC_w_sloshbox(x,flag):
    return None

dirichletConditions_ns = {0:getDBC_p_sloshbox,
                          1:getDBC_u_sloshbox,
                          2:getDBC_v_sloshbox}

def getAFBC_p_sloshbox(x,flag):
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    if not outFlowEnd and x[0] >= L[0] - 1.0e-8:
        return lambda x,t: 0.0
    if setWavesAtInflow and x[0] <= 1.0e-8:
        return waveFlux_u

def getAFBC_u_sloshbox(x,flag):
    if (x[1] == 0.0):
        return lambda x,t: 0.0
        
def getAFBC_v_sloshbox(x,flag):
    if (x[1] == 0.0):
        return lambda x,t: 0.0

def getDFBC_u_sloshbox(x,flag):
    if (x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0
    
def getDFBC_v_sloshbox(x,flag):
    if (x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0


fluxBoundaryConditions_ns = {0:'outFlow',
                             1:'outFlow',
                             2:'outFlow'}

advectiveFluxBoundaryConditions_ns =  {0:getAFBC_p_sloshbox,
                                       1:getAFBC_u_sloshbox,
                                       2:getAFBC_v_sloshbox}

diffusiveFluxBoundaryConditions_ns = {0:{},
                                      1:{1:getDFBC_u_sloshbox,2:getDFBC_u_sloshbox},
                                      2:{1:getDFBC_v_sloshbox,2:getDFBC_v_sloshbox}}



class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] >= waterLevelBase:
            return -(domain.L[1]-x[1])*rho_1*g[1]
        else:
            return -((domain.L[1]-waterLevelBase)*rho_1 +
                     (waterLevelBase-x[1])*rho_0)*g[1]

initialConditions_ns = {0:Hydrostatic_p(),
                        1:AtRest(),
                        2:AtRest()}


class Flat_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return x[1] - waterLevelBase

class Flat_H:
    def __init__(self,eps=epsFact_consrv_heaviside*dy):
        self.eps=eps
    def uOfXT(self,x,t):
        return smoothedHeaviside(self.eps,x[1] - waterLevelBase)
    
def getDBC_vof(x,flag):
    if setWavesAtInflow:
        return freeSurfaceVOF_bc


def getAFBC_vof(x,flag):
    pass
