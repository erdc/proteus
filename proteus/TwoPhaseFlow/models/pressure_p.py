from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import Pres

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
myTpFlowProblem = ct.myTpFlowProblem 
physical_parameters   = myTpFlowProblem.physical_parameters
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# DOMAIN #
domain = myTpFlowProblem.domain

# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
rho_0 = physical_parameters['densityA']
g = physical_parameters['gravity']

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
initialConditions = {0: initialConditions['pressure']}

# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: boundaryConditions['pressure_DBC']} 
advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC']}
