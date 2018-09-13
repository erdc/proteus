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
physical_parameters = myTpFlowProblem.physical_parameters
rans3p_parameters   = myTpFlowProblem.rans3p_parameters
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

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
useMetrics = rans3p_parameters['useMetrics']
epsFact_viscosity = rans3p_parameters['epsFact_viscosity']
epsFact_density = rans3p_parameters['epsFact_density']
ns_forceStrongDirichlet = rans3p_parameters['ns_forceStrongDirichlet']
weak_bc_penalty_constant = rans3p_parameters['weak_bc_penalty_constant']
useRBLES = rans3p_parameters['useRBLES']
useRANS = rans3p_parameters['useRANS']
ns_closure = rans3p_parameters['ns_closure']
useVF = rans3p_parameters['useVF']
PSTAB = rans3p_parameters['PSTAB']
USE_SUPG = rans3p_parameters['USE_SUPG']
ARTIFICIAL_VISCOSITY = rans3p_parameters['ARTIFICIAL_VISCOSITY']
cE = rans3p_parameters['cE']
cMax = rans3p_parameters['cMax']

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
CLSVOF_model=0
V_model=1
PINC_model=2
PRESSURE_model=3

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
