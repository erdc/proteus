from proteus import *
from proteus.default_p import *
from cylinder import *
from proteus.mprans import RANS3PF
from numpy import linalg as LA

name="rans3p"
LevelModelType = RANS3PF.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if useRANS >= 1:
    Closure_0_model = 5; Closure_1_model=6
    if useOnlyVF:
        Closure_0_model=2; Closure_1_model=3
    if movingDomain:
        Closure_0_model += 1; Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None
####################################################################################################
def get_repulsive_force_from_particle(Xi,Xj,Ri,Rj,force_range,stiffness_parameter):
    dij = LA.norm(Xi-Xj)
    FF = np.zeros((3,),'d')
    if dij <= Ri+Rj+force_range:
        FF = (Xi-Xj)*(Ri+Rj+force_range-dij)**(2.0)/stiffness_parameter

    return FF
def get_repulsive_force_from_left_wall(Xi,Ri,force_range,stiffness_parameter):
    Xj = np.copy(Xi)
    Xj[0] = 2*x0[0]-Xi[0]
    Rj = Ri
    return get_repulsive_force_from_particle(Xi,Xj,Ri,Rj,force_range,stiffness_parameter)
def get_repulsive_force_from_right_wall(Xi,Ri,force_range,stiffness_parameter):
    Xj = np.copy(Xi)
    Xj[0] = 2*(x0[0]+L[0])-Xi[0]
    Rj = Ri
    return get_repulsive_force_from_particle(Xi,Xj,Ri,Rj,force_range,stiffness_parameter)
def get_repulsive_force_from_bottom_wall(Xi,Ri,force_range,stiffness_parameter):
    Xj = np.copy(Xi)
    Xj[1] = 2*x0[1]-Xi[1]
    Rj = Ri
    return get_repulsive_force_from_particle(Xi,Xj,Ri,Rj,force_range,stiffness_parameter)
def get_repulsive_force_from_top_wall(Xi,Ri,force_range,stiffness_parameter):
    Xj = np.copy(Xi)
    Xj[1] = 2*(x0[1]+L[1])-Xi[1]
    Rj = Ri
    return get_repulsive_force_from_particle(Xi,Xj,Ri,Rj,force_range,stiffness_parameter)
def get_vel_from_solid_motion(C,U,omega,X):
    vel = np.copy(U)
    R = X-C
    vel[0] += -omega*R[1]
    vel[1] +=  omega*R[0]
    return vel[:2]
def distance_norm_to_disk(C,R,X):
    n = X-C
    d = LA.norm(n)
    n /= (d+1e-10)
    return d-R,n[:2]
####################################################################################################coefficient class
class My_Coefficient(RANS3PF.Coefficients):
    def __init__(self):
        RANS3PF.Coefficients.__init__(self,epsFact=epsFact_viscosity,
                                sigma=0.0,
                                rho_0 = rho_0,
                                nu_0 = nu_0,
                                rho_1 = rho_1,
                                nu_1 = nu_1,
                                g=g,
                                nd=nd,
                                ME_model=V_model,
                                PRESSURE_model=PRESSURE_model,
                                SED_model=SED_model,
                                VOF_model=VOF_model,
                                VOS_model=VOS_model,
                                LS_model=LS_model,
                                Closure_0_model=Closure_0_model,
                                Closure_1_model=Closure_1_model,
                                epsFact_density=epsFact_density,
                                stokes=False,
                                useVF=useVF,
                                useRBLES=useRBLES,
                                useMetrics=useMetrics,
                                eb_adjoint_sigma=1.0,
                                eb_penalty_constant=weak_bc_penalty_constant,
                                forceStrongDirichlet=ns_forceStrongDirichlet,
                                turbulenceClosureModel=ns_closure,
                                movingDomain=movingDomain,
                                dragAlpha=dragAlpha,
                                PSTAB=0.0,
                                nParticles=1,
                                particle_epsFact=2.0,
                                particle_alpha=1e6,
                                particle_beta=1e6,
                                particle_penalty_constant=1e16,
                                particle_sdfList=[particle_sdf_1],
                                particle_velocityList=[particle_vel_1],
                                use_sbm=USE_SBM)
        self.outputQuantDOFs=True
#     def postStep(self, t, firstStep=False):
#         RANS3PF.Coefficients.postStep(self,t,firstStep)

####################################################################################################coefficient
coefficients = My_Coefficient()
# def getDBC_u(x,flag):
#     if flag in [boundaryTags['left'],boundaryTags['right'],boundaryTags['bottom']]:
#         return lambda x,t: 0.0
# 
# def getDBC_v(x,flag):
#     if flag in [boundaryTags['left'],boundaryTags['right'],boundaryTags['bottom']]:
#         return lambda x,t: 0.0
# 
# dirichletConditions = {0:getDBC_u,
#                        1:getDBC_v}
# 
# def getAFBC_u(x,flag):
#     if flag in [boundaryTags['top']]:
#         return lambda x,t: 0.0
#     else:
#         return None
# 
# def getAFBC_v(x,flag):
#     if flag in [boundaryTags['top']]:
#         return lambda x,t: 0.0
#     else:
#         return None
# 
# def getDFBC_u(x,flag):
#     if flag in [boundaryTags['top']]:
#         return lambda x,t: 0.0
#     else:
#         return None
# 
# def getDFBC_v(x,flag):
#     if flag in [boundaryTags['top']]:
#         return lambda x,t: 0.0
#     else:
#         return None
# 
# advectiveFluxBoundaryConditions =  {0:getAFBC_u,
#                                     1:getAFBC_v}
# 
# diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
#                                    1:{1:getDFBC_v}}
####################################################################################################boundary condition
def getDBC_u(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right'],boundaryTags['top'],]:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right'],boundaryTags['top'],]:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

def getAFBC_u(x,flag):
    if flag in [boundaryTags['bottom']]:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag in [boundaryTags['bottom']]:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if flag in [boundaryTags['bottom']]:
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if flag in [boundaryTags['bottom']]:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v}}

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest()}
