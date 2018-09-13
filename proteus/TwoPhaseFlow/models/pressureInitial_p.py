from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import PresInit

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
myTpFlowProblem = ct.myTpFlowProblem 
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
PINIT_model=4
V_model=1
PRESSURE_model=3

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
#LevelModelType = PresInit.LevelModel
coefficients=PresInit.Coefficients(nd=nd,
                                   modelIndex=PINIT_model,
                                   fluidModelIndex=V_model,
                                   pressureModelIndex=PRESSURE_model)
name = "pressureInitial"

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: initialConditions['pressure']}

# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: boundaryConditions['pressure_DBC']}
advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC']}
diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['pressure_increment_DFBC']}}
