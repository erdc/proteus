from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import Pres

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
gravity ={'gx': [-9.8, 0.0, 0.0],
          'gy': [0.0, -9.8, 0.0],
          'gz': [0.0, 0.0, -9.8]}
gravity_direction = ct.physical_parameters['gravity']
g = gravity[gravity_direction]

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
PRESSURE_model=3
V_model=1
PINC_model=2

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = Pres.LevelModel
coefficients=Pres.Coefficients(modelIndex=PRESSURE_model,
                               fluidModelIndex=V_model,
                               pressureIncrementModelIndex=PINC_model,
                               useRotationalForm=False)
name = "pressure"

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: ct.pressure_init_cond()}

# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: ct.pressure_DBC} 
advectiveFluxBoundaryConditions = {0: ct.pressure_AFBC}
