from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import RANS2P
from proteus import Context

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
rho_0 = ct.physical_parameters['densityA']
nu_0 = ct.physical_parameters['viscosityA']
rho_1 = ct.physical_parameters['densityB']
nu_1 = ct.physical_parameters['viscosityB']
sigma_01 = ct.physical_parameters['surf_tension_coeff']
g = ct.physical_parameters['gravity']

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
useMetrics = ct.rans2p_parameters['useMetrics']
epsFact_viscosity = ct.rans2p_parameters['epsFact_viscosity']
epsFact_density = ct.rans2p_parameters['epsFact_density']
ns_forceStrongDirichlet = ct.rans2p_parameters['ns_forceStrongDirichlet']
weak_bc_penalty_constant = ct.rans2p_parameters['weak_bc_penalty_constant']
useRBLES = ct.rans2p_parameters['useRBLES']
useRANS = ct.rans2p_parameters['useRANS']
ns_closure = ct.rans2p_parameters['ns_closure']
useVF = ct.rans2p_parameters['useVF']

# *************************************** #
# ********** TURBULENCE MODELS ********** #
# *************************************** #
Closure_0_model = None
Closure_1_model = None

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
CLSVOF_model=1
VF_model=None
LS_model=None

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=sigma_01,
                                   rho_0=rho_0,
                                   nu_0=nu_0,
                                   rho_1=rho_1,
                                   nu_1=nu_1,
                                   g=g,
                                   nd=nd,
                                   ME_model=0,
                                   CLSVOF_model=CLSVOF_model,
                                   VF_model=VF_model,
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
                                   movingDomain=movingDomain)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
if nd==2:
    initialConditions = {0:ct.pressure_init_cond(),
                         1:ct.vel_u_init_cond(),
                         2:ct.vel_v_init_cond()}
else:
    initialConditions = {0:ct.pressure_init_cond(),
                         1:ct.vel_u_init_cond(),
                         2:ct.vel_v_init_cond(),
                         3:ct.vel_w_init_cond()}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
if nd==2:
    dirichletConditions = {0: ct.pressure_DBC,
                           1: ct.vel_u_DBC,
                           2: ct.vel_v_DBC}
    advectiveFluxBoundaryConditions = {0: ct.pressure_AFBC,
                                       1: ct.vel_u_AFBC,
                                       2: ct.vel_v_AFBC}
    diffusiveFluxBoundaryConditions = {0: {},
                                       1: {1:ct.vel_u_DFBC},
                                       2: {2:ct.vel_v_DFBC}}
else:
    dirichletConditions = {0: ct.pressure_DBC,
                           1: ct.vel_u_DBC,
                           2: ct.vel_v_DBC,
                           3: ct.vel_w_DBC}
    advectiveFluxBoundaryConditions = {0: ct.pressure_AFBC,
                                       1: ct.vel_u_AFBC,
                                       2: ct.vel_v_AFBC,
                                       3: ct.vel_w_AFBC}
    diffusiveFluxBoundaryConditions = {0: {},
                                       1: {1:ct.vel_u_DFBC},
                                       2: {2:ct.vel_v_DFBC},
                                       3: {3:ct.vel_w_DFBC}}    

