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
params = myTpFlowProblem.Parameters
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# DOMAIN #
domain = myTpFlowProblem.domain

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
rho_0 = params.physical['densityA']
nu_0 = params.physical['viscosityA']
rho_1 = params.physical['densityB']
nu_1 = params.physical['viscosityB']
sigma_01 = params.physical['surf_tension_coeff']
g = params.physical['gravity']

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
useMetrics = params.rans3p['useMetrics']
epsFact_viscosity = params.rans3p['epsFact_viscosity']
epsFact_density = params.rans3p['epsFact_density']
ns_forceStrongDirichlet = params.rans3p['ns_forceStrongDirichlet']
weak_bc_penalty_constant = params.rans3p['weak_bc_penalty_constant']
useRBLES = params.rans3p['useRBLES']
useRANS = params.rans3p['useRANS']
ns_closure = params.rans3p['ns_closure']
useVF = params.rans3p['useVF']
PSTAB = params.rans3p['PSTAB']
USE_SUPG = params.rans3p['USE_SUPG']
ARTIFICIAL_VISCOSITY = params.rans3p['ARTIFICIAL_VISCOSITY']
cE = params.rans3p['cE']
cMax = params.rans3p['cMax']

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
CLSVOF_model = params.clsvof['index']
V_model = params.rans3p['index']
PINC_model = params.pressureIncrement['index']
PRESSURE_model = params.pressure['index']

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
                                    cE=cE, cMax=cMax)

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
if nd==2:
    dirichletConditions = {0: boundaryConditions['vel_u_DBC'],
                           1: boundaryConditions['vel_v_DBC']}
    advectiveFluxBoundaryConditions =  {0: boundaryConditions['vel_u_AFBC'],
                                        1: boundaryConditions['vel_v_AFBC']}
    diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['vel_u_DFBC']},
                                       1:{1: boundaryConditions['vel_v_DFBC']}}
else:
    dirichletConditions = {0: boundaryConditions['vel_u_DBC'],
                           1: boundaryConditions['vel_v_DBC'],
                           2: boundaryConditions['vel_w_DBC']}
    advectiveFluxBoundaryConditions =  {0: boundaryConditions['vel_u_AFBC'],
                                        1: boundaryConditions['vel_v_AFBC'],
                                        2: boundaryConditions['vel_w_AFBC']}
    diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['vel_u_DFBC']},
                                       1:{1: boundaryConditions['vel_v_DFBC']},
                                       2:{2: boundaryConditions['vel_w_DFBC']}}
