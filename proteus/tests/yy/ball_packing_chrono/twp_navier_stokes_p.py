from proteus import *
from proteus.default_p import *
from rockyRiver import *
from proteus.mprans import RANS3PF

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


mylist_of_sdfs=[];
# for i in range(numPars):
#     mylist_of_sdfs.append(particle_sdf)

# coefficients = RANS3PF.Coefficients(epsFact=epsFact_viscosity,
#                                     sigma=0.0,
#                                     rho_0 = rho_0,
#                                     nu_0 = nu_0,
#                                     rho_1 = rho_1,
#                                     nu_1 = nu_1,
#                                     g=g,
#                                     nd=nd,
#                                     ME_model=V_model,
#                                     PRESSURE_model=PRESSURE_model,
#                                     SED_model=SED_model,
#                                     VOF_model=VOF_model,
#                                     VOS_model=VOS_model,
#                                     LS_model=LS_model,
#                                     Closure_0_model=Closure_0_model,
#                                     Closure_1_model=Closure_1_model,
# #                                     epsFact_density=epsFact_density,
#                                     stokes=False,
#                                     useVF=useVF,
#                                     useRBLES=useRBLES,
#                                     useMetrics=useMetrics,
#                                     eb_adjoint_sigma=1.0,
#                                     eb_penalty_constant=weak_bc_penalty_constant,
#                                     forceStrongDirichlet=ns_forceStrongDirichlet,
#                                     turbulenceClosureModel=ns_closure,
#                                     movingDomain=movingDomain,
#                                     dragAlpha=dragAlpha,
#                                     PSTAB=1.0,
#                                     nParticles=numPars,
#                                     particle_epsFact=2.0,
#                                     particle_alpha=1e3,
#                                     particle_beta=1e3,
#                                     particle_penalty_constant=1e6,
#                                     particle_sdfList=[],
#                                     particle_velocityList=[],
#                                     granular_sdf_Calc=my_granular_sdf_Calc,
#                                     granular_vel_Calc= my_granular_vel_Calc,
#                                     use_sbm=USE_SBM)

class My_coefficients(RANS3PF.Coefficients):
    def __init__(self):
        RANS3PF.Coefficients.__init__(self,
                                    epsFact=epsFact_viscosity,
                                    sigma=0.0,
                                    USE_SUPG = use_supg,#####
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
#                                     epsFact_density=epsFact_density,
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
                                    PSTAB=1.0,
                                    nParticles=numPars,
                                    particle_epsFact=1.0,
                                    particle_alpha=1e3,
                                    particle_beta=1e3,
                                    particle_penalty_constant=1e6,
                                    particle_sdfList=[],
                                    particle_velocityList=[],
                                    use_sbm=USE_SBM,
                                    use_ball_as_particle = use_ball_as_particle,
                                    ball_center=ball_center,
                                    ball_radius=ball_radius,
                                    ball_velocity=ball_velocity,
                                    ball_angular_velocity=ball_angular_velocity,
                                    )

#         self.ARTIFICIAL_VISCOSITY = 2#0: no art viscosity, 1: shock capturing (default), 2: entropy viscosity
    def attachModels(self, modelList):
        RANS3PF.Coefficients.attachModels(self,modelList)
        
        self.model=modelList[0]##momentum module
        self.model.q['q_phi_solid'] = self.q_phi_solid
    
    def preStep(self, t, firstStep=False):
        
        self.particle_netForces[:] = 0.0
        self.particle_netMoments[:] = 0.0

        self.model.u_dof_old_old[:] = self.model.u_dof_old[:]
        self.model.v_dof_old_old[:] = self.model.v_dof_old[:]
    
        self.model.u_dof_old[:] = self.model.u[0].dof
        self.model.v_dof_old[:] = self.model.u[1].dof
        
#         if (firstStep):
#             self.model.q[('velocityStar', 0)][:] = self.model.q[('velocity', 0)]
#         else:
#             if self.model.timeIntegration.timeOrder == 1:
#                 r = 1.
#             else:
#                 r = self.model.timeIntegration.dt / self.model.timeIntegration.dt_history[0]
#             self.model.q[('velocityStar', 0)][:] = (1 + r) * self.model.q[('velocity', 0)] - r * self.model.q[('velocityOld', 0)]

        self.model.q[('velocityStar', 0)][:] = self.model.q[('velocity', 0)]
        self.model.q[('velocityOld', 0)][:] = self.model.q[('velocity', 0)]

        self.model.dt_last = self.model.timeIntegration.dt

    def postStep(self, t, firstStep=False):
        ##ball center is updated in Chrono 
        RANS3PF.Coefficients.postStep(self, t, firstStep)

coefficients = My_coefficients()
#===============================================================================
# BC
#===============================================================================

def getDBC_u(x,flag):
    if flag in [ boundaryTags['left'] ,boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom']]:
        return lambda x,t: 0.0
    else:
        None

def getDBC_v(x,flag):
    if flag in [ boundaryTags['left'] ,boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom']]:
        return lambda x,t: 0.0
    else:
        None

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

def getAFBC_u(x,flag):
    if flag in [ boundaryTags['left'] ,boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom']]:
        return None
    else:
        return None#MUST be None; consistent with pressure increment module;

def getAFBC_v(x,flag):
      if flag in [ boundaryTags['left'] ,boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom']]:
          return None
      else:
          return None#MUST be None; consistent with pressure increment module;

def getDFBC_u(x,flag):
  if flag in [ boundaryTags['left'] ,boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom']]:
      return None
  else:
      return lambda x,t: 0.0

def getDFBC_v(x,flag):
  if flag in [ boundaryTags['left'] ,boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom']]:
      return None
  else:
      return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v}}

class AtRest:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest()}
