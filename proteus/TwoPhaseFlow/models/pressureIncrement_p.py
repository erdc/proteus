from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import PresInc

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
myTpFlowProblem = ct.myTpFlowProblem 
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# DOMAIN #
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
rho_0 = pparams['densityA']
rho_1 = pparams['densityB']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
PINC_model = mparams.pressureIncrement['index']
V_model = mparams.rans3p['index']

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = PresInc.LevelModel
coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
                                  rho_s_min = (1.0-1.0e-8)*rho_0,
                                  nd = nd,
                                  modelIndex=PINC_model,
                                  fluidModelIndex=V_model,
                                  fixNullSpace=False)
name = "pressureIncrement"

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
initialConditions = {0: initialConditions['pressure_increment']}

# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['pressure_increment_DBC']}
    advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_increment_AFBC']}
    diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['pressure_increment_DFBC']}}
else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].pInc_dirichlet.init_cython()}
    advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].pInc_advective.init_cython()}
    diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].pInc_diffusive.init_cython()}}
