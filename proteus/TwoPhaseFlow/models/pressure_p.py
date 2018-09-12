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
g = ct.physical_parameters['gravity']

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
