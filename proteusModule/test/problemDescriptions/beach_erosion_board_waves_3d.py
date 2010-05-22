from math import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from pyadh import Domain
from pyadh import TransportCoefficients
from pyadh.TransportCoefficients import TwophaseNavierStokes_ST_LS_SO
from pyadh import LinearAlgebraTools
from pyadh import MeshTools
import beach_erosion_board_3dDomain
#numerics options
useVANS2P=True
useNCLS  =True
useVOF   =True
useRDLS  =True
useMCorr =True

applyCorrection  =True#True#False
applyRedistancing=True#right now the solution is more symmetric without redistancing
rdtimeIntegration='osher'#'osher-psitc'#tte
sloshbox_quad_order = 3 #need 4 for quadratics

checkMass = True
useBackwardEuler=True#True#False
useBackwardEuler_ls=True#False

nonlinearSolverNorm = LinearAlgebraTools.l2NormAvg

lag_ns_subgridError=True
lag_ls_shockCapturing=True
ls_shockCapturingFactor=0.33#0.1
rd_shockCapturingFactor=0.33#0.1
ns_shockCapturingFactor=0.33

rd_freezeLS = True
usePETSc=True#False#True
nOverlap = 1
partitioningType = MeshTools.MeshParallelPartitioningTypes.element
#-P"-ksp_type bcgsl -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_ksp_type preonly -sub_pc_type lu -sub_pc_factor_mat_solver_package spooles -ksp_atol 1.0e-7 -ksp_rtol 1.0e-7 "
#
useStokes=False
VOF = False
nonconsrv = True
if  VOF:
    nonconsrv = True
spaceOrder=1

#spatial domain
nd = 3

domain = beach_erosion_board_3dDomain.beach_erosion_board_3d(domainHeightPad=0.2,
                                                             inflowLength=3.0,#0.45,3.0
                                                             beachLength=4.48,
                                                             beachSlope =0.1,
                                                             h_c=0.054,#0.054,
                                                             h_s=0.081,
                                                             s  =3.0,
                                                             B  = 0.0896,
                                                             s_2= -0.2,
                                                             backStepFraction=0.25,
                                                             outflowLength=0.5,
                                                             width = 0.3)
nnx=3; nny=3
nLevels = 1
height=domain.L[2]
length=domain.L[0]
width =domain.L[1]
L = domain.L
constraintFactor = 0.2#0.1#0.2
volumeConstraint =  (((constraintFactor*height)**3)/6.0,)
dx =0.05*height ; dy = dx ; dz = dx; he = dx;
triangleOptions="VpAq1.25ena%f" % (volumeConstraint)
#triangleOptions="VpAq1.25en" 


domain.writeAsymptote("beach_erosion_board_3d")
domain.writePoly("beach_erosion_board_3d")



#wave specification
runWaveProblem=True

#source zone for generating waves
#background water level
waterLevelBase = 0.529#domain.bathymetry([domain.x_be])+domain.h_s#backHeight#0.529#should be 0.529 #m
waveHeight=0.107#0.5#default 0.107#m
wavePeriod= 1.549   #[s]
waveCelerity  = 1.21#[m/s] don't know
waveNumber= 1.04/waterLevelBase #[1/m]
waveLength= 2.0*pi/waveNumber
waveFrequency = 2.0*pi/wavePeriod
source_height=max(5*dy,0.5*waterLevelBase)#max(5*dy,L[1])#max(5*dy,0.25*waterLevelBase)
source_zm    =domain.x[2]+0.3*waterLevelBase
#source_ym    =domain.x[1]+0.5*L[1]#waterLevelBase
source_xm    =0.75*domain.inflowLength#0.5*domain.inflowLength
source_length =max(3.0*dx,0.1)#better max(3.0*dx,0.2)
Omega_s=[[source_xm,source_xm+source_length],
         [0,width],
         [source_zm-0.5*source_height,source_zm+0.5*source_height]]
    
sourceVolume = (Omega_s[0][1]-Omega_s[0][0])*(Omega_s[1][1]-Omega_s[1][0])*(Omega_s[2][1]-Omega_s[2][0])
waveFlag= 0 #0 -- monochromatic
            #1 -- second order stokes
            #2 -- solitary wave
            #-1 -- no internal wave
setWavesAtInflow = False#True
if setWavesAtInflow:
    waveFlag = -1

topOpen = True#False


#boundary conditions
#slip on bottom, 
#outflow on sides
#top open with grate --> p = 0, u = 0, v is free
 

#
#numerical fudge factors
epsFact = 1.5
epsFact_density = 1.5#0.33#1.5#1.5
epsFact_viscosity = 1.5#0.33#1.5
epsFact_redistance = 0.33#1.5
epsFact_curvature=0.001
epsFact_source = 0.5
epsFact_vof=1.5
epsFact_consrv_heaviside = 1.5
epsFact_consrv_dirac = 1.5
epsFact_consrv_diffusion = 10.0#1.5

#physical parameters
#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5
#gravity
g=[0.0,0.0,-9.8]
#mwf play with density and viscosity ratios
#rho_1 = rho_0/20.
#nu_1  = nu_0*0.5
turbulenceClosureFlag = 1#Smagorinsky
smagorinskyConstant_0 = 0.1
smagorinskyConstant_1 = 0.5
#surface tension
sigma_01 = 0.0#72.8e-3

import math
T =0.3#wavePeriod*5.0#*10.0#10.#*10.0#2.5#3.0#10.0
nDTout = 3#math.ceil(T/wavePeriod)*10+1#51#101#None#2#None, can't be 1
runCFL = 0.33

dt_init = min(0.01*wavePeriod,T)
#set up absorbing boundary at left side using volume averaged ns formulation
spongeLayerWidth = max(2.0*dx,0.3)#0.2*domain.inflowLength)
useSpongeLayer = True#True#False
spongeGrainSize= 0.01
spongePorosity = 0.5
killNonlinearDragInSpongeLayer = True#True
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
bcsTimeDependent = setWavesAtInflow 
def getDBC_p_ns(x,flag):
    if flag == domain.boundaryTags['top']:
        return lambda x,t: 0.0

def getDBC_u_ns(x,flag):
    if flag  == domain.boundaryTags['top']:
        return lambda x,t: 0.0

    
def getDBC_v_ns(x,flag):
    if flag  == domain.boundaryTags['top']:
        return lambda x,t: 0.0


def getDBC_w_ns(x,flag):
    if not topOpen and flag  == domain.boundaryTags['top']:
        return lambda x,t: 0.0

dirichletConditions_ns = {0:getDBC_p_ns,
                          1:getDBC_u_ns,
                          2:getDBC_v_ns,
                          3:getDBC_w_ns}

slip_walls = [domain.boundaryTags['bottom'],domain.boundaryTags['front'],domain.boundaryTags['back']]
outflow_walls = [domain.boundaryTags['left'],domain.boundaryTags['right']]
def getAFBC_p_ns(x,flag):
    if flag in slip_walls:
        return lambda x,t: 0.0
    if not topOpen and flag  == domain.boundaryTags['top']:
        return lambda x,t: 0.0
    
def getAFBC_u_ns(x,flag):
    if flag in slip_walls:
        return lambda x,t: 0.0
        
def getAFBC_v_ns(x,flag):
    if flag in slip_walls:
        return lambda x,t: 0.0

def getAFBC_w_ns(x,flag):
    if flag in slip_walls:
        return lambda x,t: 0.0

def getDFBC_u_ns(x,flag):
    #go ahead and enforce no diffusive flux on outflow too for VANS2P
    #need parallel boundaries as well
    return lambda x,t: 0.0
    if flag in slip_walls:# or flag in [domain.boundaryTags['top']]:
        return lambda x,t: 0.0
    if flag in outflow_walls:
        return lambda x,t: 0.0
def getDFBC_v_ns(x,flag):
    #need parallel boundaries as well
    return lambda x,t: 0.0
    if flag in  slip_walls:# or flag in [domain.boundaryTags['top']]:
        return lambda x,t: 0.0
    #go ahead and enforce no diffusive flux on outflow too for VANS2P
    if flag in outflow_walls:
        return lambda x,t: 0.0
def getDFBC_w_ns(x,flag):
    #need parallel boundaries as well
    return lambda x,t: 0.0
    if flag in  slip_walls:# or flag in [domain.boundaryTags['top']]:
        return lambda x,t: 0.0
    #go ahead and enforce no diffusive flux on outflow too for VANS2P
    if flag in outflow_walls:
        return lambda x,t: 0.0

fluxBoundaryConditions_ns = {0:'outFlow',
                             1:'outFlow',
                             2:'outFlow',
                             3:'outFlow'}

advectiveFluxBoundaryConditions_ns =  {0:getAFBC_p_ns,
                                       1:getAFBC_u_ns,
                                       2:getAFBC_v_ns,
                                       3:getAFBC_w_ns}

diffusiveFluxBoundaryConditions_ns = {0:{},
                                      1:{1:getDFBC_u_ns},#,2:getDFBC_u_ns,3:getDFBC_u_ns},
                                      2:{2:getDFBC_v_ns},#3:getDFBC_v_ns},
                                      3:{3:getDFBC_w_ns}}


#try to vary free surface height at boundary

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[2] >= waterLevelBase:
            return -(domain.L[2]-x[2])*rho_1*g[2]
        else:
            return -((domain.L[2]-waterLevelBase)*rho_1 +
                     (waterLevelBase-x[2])*rho_0)*g[2]

initialConditions_ns = {0:Hydrostatic_p(),
                        1:AtRest(),
                        2:AtRest(),
                        3:AtRest()}


class Flat_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return x[2] - waterLevelBase

class Flat_H:
    def __init__(self,eps=epsFact_consrv_heaviside*dy):
        self.eps=eps
    def uOfXT(self,x,t):
        return smoothedHeaviside(self.eps,x[2] - waterLevelBase)
    
def getDBC_vof(x,flag):
    if setWavesAtInflow:
        return freeSurfaceVOF_bc
    if flag in [domain.boundaryTags['top']]:
        return lambda x,t: 1.0
    #enforce air on right outflow boundary?
    if flag in [domain.boundaryTags['right']]:
        return lambda x,t: 1.0
    #
    if flag in [domain.boundaryTags['left']]:
        return lambda x,t: smoothedHeaviside(epsFact*dz,x[2]-waterLevelBase)
def getAFBC_vof(x,flag):
    #pass
    if flag == 0: #internal boundary
        return lambda x,t: 0.0
   
    if flag in  slip_walls:# or flag in [domain.boundaryTags['top']]:
        return lambda x,t: 0.0
