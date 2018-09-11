from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import PresInc

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
rho_0 = ct.physical_parameters['densityA']
rho_1 = ct.physical_parameters['densityB']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
PINC_model=2
V_model=1

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
name = "pressureincrement"

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
initialConditions = {0: ct.pressure_increment_init_cond()}

# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
dirichletConditions = {0: ct.pressure_increment_DBC}
advectiveFluxBoundaryConditions = {0: ct.pressure_increment_AFBC}
diffusiveFluxBoundaryConditions = {0:{0: ct.pressure_increment_DFBC}}
