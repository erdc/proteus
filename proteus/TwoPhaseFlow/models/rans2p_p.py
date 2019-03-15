from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import RANS2P
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
<<<<<<< HEAD
myTpFlowProblem = ct.myTpFlowProblem 
physical_parameters = myTpFlowProblem.physical_parameters
rans2p_parameters   = myTpFlowProblem.rans2p_parameters
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
=======
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

myTpFlowProblem = ct.myTpFlowProblem 
IC = myTpFlowProblem.initialConditions
boundaryConditions = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain
>>>>>>> TwoPhaseFlow

# DOMAIN #
domain = myTpFlowProblem.domain

<<<<<<< HEAD
# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
rho_0 = physical_parameters['densityA']
nu_0 = physical_parameters['viscosityA']
rho_1 = physical_parameters['densityB']
nu_1 = physical_parameters['viscosityB']
sigma_01 = physical_parameters['surf_tension_coeff']
g = physical_parameters['gravity']

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
useMetrics = rans2p_parameters['useMetrics']
epsFact_viscosity = rans2p_parameters['epsFact_viscosity']
epsFact_density = rans2p_parameters['epsFact_density']
ns_forceStrongDirichlet = rans2p_parameters['ns_forceStrongDirichlet']
weak_bc_penalty_constant = rans2p_parameters['weak_bc_penalty_constant']
useRBLES = rans2p_parameters['useRBLES']
useRANS = rans2p_parameters['useRANS']
ns_closure = rans2p_parameters['ns_closure']
useVF = rans2p_parameters['useVF']
=======
params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
myparams = params.Models.rans2p # model parameters
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
>>>>>>> TwoPhaseFlow

# *************************************** #
# ********** TURBULENCE MODELS ********** #
# *************************************** #
Closure_0_model = None
Closure_1_model = None

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
<<<<<<< HEAD
CLSVOF_model=1
VF_model=None
LS_model=None

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
=======
ME_model = mparams.rans2p.index
assert ME_model is not None, 'rans2p model index was not set!'
CLSVOF_model = mparams.clsvof.index
VF_model = mparams.vof.index
LS_model = mparams.ncls.index

# for absorption zones (defined as regions)
if hasattr(domain, 'porosityTypes'):
    porosityTypes = domain.porosityTypes
    dragAlphaTypes = domain.dragAlphaTypes
    dragBetaTypes = domain.dragBetaTypes
    epsFact_solid = domain.epsFact_solid
else:
    porosityTypes = None
    dragAlphaTypes = None
    dragBetaTypes = None
    epsFact_solid = None

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=myparams.epsFact_viscosity,
>>>>>>> TwoPhaseFlow
                                   sigma=sigma_01,
                                   rho_0=rho_0,
                                   nu_0=nu_0,
                                   rho_1=rho_1,
                                   nu_1=nu_1,
                                   g=g,
                                   nd=nd,
<<<<<<< HEAD
                                   ME_model=0,
=======
                                   ME_model=ME_model,
>>>>>>> TwoPhaseFlow
                                   CLSVOF_model=CLSVOF_model,
                                   VF_model=VF_model,
                                   LS_model=LS_model,
                                   Closure_0_model=Closure_0_model,
                                   Closure_1_model=Closure_1_model,
<<<<<<< HEAD
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
=======
                                   epsFact_density=myparams.epsFact_density,
                                   stokes=myparams.stokes,
                                   useVF=myparams.useVF,
                                   useRBLES=myparams.useRBLES,
                                   useMetrics=myparams.useMetrics,
                                   eb_adjoint_sigma=myparams.eb_adjoint_sigma,
                                   eb_penalty_constant=myparams.weak_bc_penalty_constant,
                                   forceStrongDirichlet=myparams.ns_forceStrongDirichlet,
                                   turbulenceClosureModel=myparams.ns_closure,
                                   movingDomain=movingDomain,
                                   porosityTypes=porosityTypes,
                                   dragAlphaTypes=dragAlphaTypes,
                                   dragBetaTypes=dragBetaTypes,
                                   epsFact_solid=epsFact_solid,
                                   barycenters=domain.barycenters
)
>>>>>>> TwoPhaseFlow

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
<<<<<<< HEAD
if nd==2:
    initialConditions = {0: initialConditions['pressure'],
                         1: initialConditions['vel_u'],
                         2: initialConditions['vel_v']}
else:
    initialConditions = {0: initialConditions['pressure'],
                         1: initialConditions['vel_u'],
                         2: initialConditions['vel_v'],
                         3: initialConditions['vel_w']}
=======
initialConditions = {0: IC['pressure'],
                     1: IC['vel_u'],
                     2: IC['vel_v']}
if nd == 3:
    initialConditions[3] = IC['vel_w']
>>>>>>> TwoPhaseFlow

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
<<<<<<< HEAD
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
    
=======

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
        advectiveFluxBoundaryConditions[3] = lambda x, flag: domain.bc[flag].w_advective.init_cython()
        diffusiveFluxBoundaryConditions[3] = {3: lambda x, flag: domain.bc[flag].w_diffusive.init_cython()}
>>>>>>> TwoPhaseFlow
