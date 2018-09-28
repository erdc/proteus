from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import RANS2P
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

myTpFlowProblem = ct.myTpFlowProblem 
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain

# DOMAIN #
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
rho_0 = pparams['densityA']
nu_0 = pparams['viscosityA']
rho_1 = pparams['densityB']
nu_1 = pparams['viscosityB']
sigma_01 = pparams['surf_tension_coeff']
g = pparams['gravity']

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
useMetrics = mparams.rans2p['useMetrics']
epsFact_viscosity = mparams.rans2p['epsFact_viscosity']
epsFact_density = mparams.rans2p['epsFact_density']
ns_forceStrongDirichlet = mparams.rans2p['ns_forceStrongDirichlet']
weak_bc_penalty_constant = mparams.rans2p['weak_bc_penalty_constant']
useRBLES = mparams.rans2p['useRBLES']
useRANS = mparams.rans2p['useRANS']
ns_closure = mparams.rans2p['ns_closure']
useVF = mparams.rans2p['useVF']

# *************************************** #
# ********** TURBULENCE MODELS ********** #
# *************************************** #
Closure_0_model = None
Closure_1_model = None

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = mparams.rans2p['index']
assert ME_model is not None, 'rans2p model index was not set!'
CLSVOF_model = mparams.clsvof['index']
VF_model = mparams.vof['index']
LS_model = mparams.ncls['index']

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   sigma=sigma_01,
                                   rho_0=rho_0,
                                   nu_0=nu_0,
                                   rho_1=rho_1,
                                   nu_1=nu_1,
                                   g=g,
                                   nd=nd,
                                   ME_model=ME_model,
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
                                   eb_adjoint_sigma=1.,
                                   eb_penalty_constant=weak_bc_penalty_constant,
                                   forceStrongDirichlet=ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ns_closure,
                                   movingDomain=movingDomain)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: initialConditions['pressure'],
                     1: initialConditions['vel_u'],
                     2: initialConditions['vel_v']}
if nd == 3:
    initialConditions[3] = initialConditions['vel_w']

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #

if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    if nd==2:
        dirichletConditions = {0: boundaryConditions['pressure_DBC'],
                               1: boundaryConditions['vel_u_DBC'],
                               2: boundaryConditions['vel_v_DBC']}
        advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC'],
                                           1: boundaryConditions['vel_u_AFBC'],
                                           2: boundaryConditions['vel_v_AFBC']}
        diffusiveFluxBoundaryConditions = {0: {},
                                           1: {1: boundaryConditions['vel_u_DFBC']},
                                           2: {2: boundaryConditions['vel_v_DFBC']}}
    else:
        dirichletConditions = {0: boundaryConditions['pressure_DBC'],
                               1: boundaryConditions['vel_u_DBC'],
                               2: boundaryConditions['vel_v_DBC'],
                               3: boundaryConditions['vel_w_DBC']}
        advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC'],
                                           1: boundaryConditions['vel_u_AFBC'],
                                           2: boundaryConditions['vel_v_AFBC'],
                                           3: boundaryConditions['vel_w_AFBC']}
        diffusiveFluxBoundaryConditions = {0: {},
                                           1: {1: boundaryConditions['vel_u_DFBC']},
                                           2: {2: boundaryConditions['vel_v_DFBC']},
                                           3: {3: boundaryConditions['vel_w_DFBC']}}
else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
                           1: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                           2: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}
    advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython(),
                                       1: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                       2: lambda x, flag: domain.bc[flag].v_advective.init_cython()}
    diffusiveFluxBoundaryConditions = {0: {},
                                       1: {1:lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                       2: {2:lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}
    if nd == 3:
        dirichletConditions[3] = lambda x, flag: domain.bc[flag].w_dirichlet.init_cython()
        advectiveFLuxBoundaryConditions[3] = lambda x, flag: domain.bc[flag].w_advective.init_cython()
        diffusiveFluxBoundaryConditions[3] = {3: lambda x, flag: domain.bc[flag].w_diffusive.init_cython()}
