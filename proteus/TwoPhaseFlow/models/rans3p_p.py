from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import RANS3PF

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain
forceTerms = None
MULTIPLY_EXTERNAL_FORCE_BY_DENSITY=0

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
nu_0 = pparams['kinematicViscosityA']
rho_1 = pparams['densityB']
nu_1 = pparams['kinematicViscosityB']
sigma_01 = pparams['surf_tension_coeff']
g = pparams['gravity']

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
if mparams.rans3p.forceTerms is not None:
    forceTerms = mparams.rans3p.forceTerms
    MULTIPLY_EXTERNAL_FORCE_BY_DENSITY=1
useMetrics = mparams.rans3p['useMetrics']
epsFact_viscosity = mparams.rans3p['epsFact_viscosity']
epsFact_density = mparams.rans3p['epsFact_density']
ns_forceStrongDirichlet = mparams.rans3p['ns_forceStrongDirichlet']
weak_bc_penalty_constant = mparams.rans3p['weak_bc_penalty_constant']
useRBLES = mparams.rans3p['useRBLES']
useRANS = mparams.rans3p['useRANS']
ns_closure = mparams.rans3p['ns_closure']
useVF = mparams.rans3p['useVF']
PSTAB = mparams.rans3p['PSTAB']
USE_SUPG = mparams.rans3p['USE_SUPG']
ARTIFICIAL_VISCOSITY = mparams.rans3p['ARTIFICIAL_VISCOSITY']
INT_BY_PARTS_PRESSURE = mparams.rans3p['INT_BY_PARTS_PRESSURE']
cE = mparams.rans3p['cE']
cMax = mparams.rans3p['cMax']

# *************************************** #
# ********** TURBULENCE MODELS ********** #
# *************************************** #
Closure_0_model = None
Closure_1_model = None

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
VOF_model=None
LS_model=None
RD_model=None
MCORR_model=None
SED_model=None
VOS_model=None
CLSVOF_model = mparams.clsvof['index']
V_model = mparams.rans3p['index']
PINC_model = mparams.pressureIncrement['index']
PRESSURE_model = mparams.pressure['index']

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = RANS3PF.LevelModel
coefficients = RANS3PF.Coefficients(epsFact=epsFact_viscosity,
                                    sigma=sigma_01,
                                    rho_0 = rho_0,
                                    nu_0 = nu_0,
                                    rho_1 = rho_1,
                                    nu_1 = nu_1,
                                    g=g,
                                    nd=nd,
                                    ME_model=V_model,
                                    PRESSURE_model=PRESSURE_model,
                                    SED_model=SED_model,
                                    CLSVOF_model=CLSVOF_model,                                    
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
                                    PSTAB=PSTAB,
                                    USE_SUPG=USE_SUPG,
                                    ARTIFICIAL_VISCOSITY=ARTIFICIAL_VISCOSITY,
                                    INT_BY_PARTS_PRESSURE=INT_BY_PARTS_PRESSURE,
                                    cE=cE, cMax=cMax,
                                    MULTIPLY_EXTERNAL_FORCE_BY_DENSITY=MULTIPLY_EXTERNAL_FORCE_BY_DENSITY)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
if nd==2:
    initialConditions = {0: initialConditions['vel_u'],
                         1: initialConditions['vel_v']}
else:
    initialConditions = {0: initialConditions['vel_u'],
                         1: initialConditions['vel_v'],
                         2: initialConditions['vel_w']}

# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #

if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['vel_u_DBC'],
                           1: boundaryConditions['vel_v_DBC']}
    advectiveFluxBoundaryConditions = {0: boundaryConditions['vel_u_AFBC'],
                                       1: boundaryConditions['vel_v_AFBC']}
    diffusiveFluxBoundaryConditions = {0: {0: boundaryConditions['vel_u_DFBC']},
                                       1: {1: boundaryConditions['vel_v_DFBC']}}
    if nd == 3:
        dirichletConditions[2] = boundaryConditions['vel_w_DBC']
        advectiveFluxBoundaryConditions[2] = boundaryConditions['vel_w_AFBC']
        diffusiveFluxBoundaryConditions[2] = {2: boundaryConditions['vel_w_DFBC']}
else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                           1: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}
    advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                       1: lambda x, flag: domain.bc[flag].v_advective.init_cython()}
    diffusiveFluxBoundaryConditions = {0: {0:lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                       1: {1:lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}
    if nd == 3:
        dirichletConditions[2] = lambda x, flag: domain.bc[flag].w_dirichlet.init_cython()
        advectiveFLuxBoundaryConditions[2] = lambda x, flag: domain.bc[flag].w_advective.init_cython()
        diffusiveFluxBoundaryConditions[2] = {2: lambda x, flag: domain.bc[flag].w_diffusive.init_cython()}
