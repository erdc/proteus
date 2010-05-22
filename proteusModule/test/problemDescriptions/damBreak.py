from math import *
import pyadh.MeshTools
from boxCylinderDomain2D import *

applyCorrection=True
applyRedistancing=True
#rdtimeIntegration='newton'
#freezeLevelSet=False
rdtimeIntegration='osher'
freezeLevelSet=True#False
damBreak_quad_order = 3
nd = 2
useBackwardEuler=True
useBackwardEuler_ls=True

dt_init=1.0e-4
T=20.0
nDTout=100
runCFL = 0.33

lag_ns_subgridError=True
lag_ls_shockCapturing=True

ns_shockCapturingFactor=0.9
ls_shockCapturingFactor=0.9
vof_shockCapturingFactor=0.9
rd_shockCapturingFactor=0.9

#epsilons for Heaviside/Dirac/etc smoothing
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.75
epsFact_curvature=1.5
epsFact_consrv_heaviside=1.5
epsFact_consrv_dirac=1.5
epsFact_consrv_diffusion=10.0
epsFact_vof=1.5

usePETSc=False

spaceOrder=1

L = (10.0,2.0,1.0)
domain = boxCylinder2D(height=L[1],length=L[0],
                       radius=0.25*L[1],cross_section=circular_cross_section,center=[0.4*L[0],0.5*L[1]],
                       boxStart=0.2*L[0],boxWidth=0.05*L[0],boxHeight=0.25*L[1])
domain.writeAsymptote("boxCylinder2D")
domain.writePoly("boxCylinder2D")

useStokes=False

he = L[1]/10.0

triangleOptions="pAq30Dena%f" % (0.5*(he**2),)
nLevels = 1

useShock=False#True
shock_x = 0.5*L[0]
shock_y = 0.5*L[1]

waterLevel = 0.5*L[1]
slopeAngle = 0.0#0.5*(pi/2.0)

#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5
sigma_01=0.0#72.8e-3
#gravity
g=[0.0,-9.8,0.0]
RE = 1.0e6
U = nu_0*RE/(2.0*0.25*L[1])

walls = [domain.boundaryTags['right'],
         domain.boundaryTags['bottom'],
         domain.boundaryTags['left'],
         domain.boundaryTags['obstacle']]

def getDBC_p_damBreak(x,flag):
    if flag == domain.boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == domain.boundaryTags['right']:
        if x[1] <= waterLevel:
            return lambda x,t: -(
                (L[1]-waterLevel)*rho_1 + 
                (waterLevel-x[1])*rho_0
                )*g[1]
        else:
            return lambda x,t: -(L[1]-x[1])*rho_1*g[1]
#        return lambda x,t: -(L[1]-x[1])*g[1]*rho_1

def getDBC_u_damBreak(x,flag):
    if flag == domain.boundaryTags['left']:
        print "left ",x[0]
        return lambda x,t: U
    if flag == domain.boundaryTags['right']:
        print "right",x[0],
        return None
    if flag in walls:
        return lambda x,t: 0.0

def getDBC_v_damBreak(x,flag):
    if flag == domain.boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == domain.boundaryTags['right']:
        return None
    if flag in walls:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_damBreak,
                       1:getDBC_u_damBreak,
                       2:getDBC_v_damBreak}

def getAFBC_p_damBreak(x,flag):
    if flag == domain.boundaryTags['left']:
        return lambda x,t: -U
    if flag == domain.boundaryTags['right']:
        return None
    if flag in walls:
        return lambda x,t: 0.0

def getAFBC_u_damBreak(x,flag):
    if flag == domain.boundaryTags['right']:
        return None
    return None

def getAFBC_v_damBreak(x,flag):
    if flag == domain.boundaryTags['right']:
        return None
    return None

def getDFBC_u_damBreak(x,flag):
    if flag == domain.boundaryTags['top']:
        return lambda x,t: 0.0
    
def getDFBC_v_damBreak(x,flag):
    if flag == domain.boundaryTags['top']:
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_damBreak,
                                    1:getAFBC_u_damBreak,
                                    2:getAFBC_v_damBreak}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_damBreak,2:getDFBC_u_damBreak},
                                   2:{1:getDFBC_v_damBreak,2:getDFBC_v_damBreak}}

def shockSignedDistance(x):
    phi_x = shock_x - x[0]
    phi_y = x[1] - shock_y
    phi2 = sqrt((x[0]-shock_x)**2 + (x[1] - shock_y)**2)
    if x[0] < shock_x:
        if x[1] < shock_y:
            return phi_x
        else:
            return phi2
    else:
        if x[1] > shock_y:
            return phi_y
        elif fabs(phi_y) < fabs(phi_x):
            return phi_y
        else:
            return phi_x

class PerturbedSurface_p:
    def __init__(self,waterLevel,slopeAngle,shock_x,shock_y):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
        self.shock_x=shock_x
        self.shock_y=shock_y
    def uOfXT(self,x,t):
        if useShock:
            if x[0] < shock_x:
                return -rho_1*g[1]*(L[1] - x[1])
            else:
                if x[1] > shock_y:
                    return -rho_1*g[1]*(L[1] - x[1])
                else:
                    return -rho_1*g[1]*(L[1] - shock_y)-rho_0*g[1]*(shock_y - x[1])
        else:
            if x[1] > waterLevel:
                return -rho_1*g[1]*(L[1] - x[1])
            else:
                return -rho_1*g[1]*(L[1] - waterLevel)-rho_0*g[1]*(waterLevel - x[1])
            

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:PerturbedSurface_p(waterLevel,slopeAngle,shock_x,shock_y),
                     1:AtRest(),
                     2:AtRest()}

restrictFineSolutionToAllMeshes=False
parallelPartitioningType = pyadh.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 2
