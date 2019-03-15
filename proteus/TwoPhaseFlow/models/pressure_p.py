from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import Pres

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
<<<<<<< HEAD
myTpFlowProblem = ct.myTpFlowProblem 
physical_parameters   = myTpFlowProblem.physical_parameters
=======

myTpFlowProblem = ct.myTpFlowProblem 
>>>>>>> TwoPhaseFlow
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# DOMAIN #
domain = myTpFlowProblem.domain

<<<<<<< HEAD
# ***************************************** #
# ********** PHYSICAL PARAMETERS ********** #
# ***************************************** #
rho_0 = physical_parameters['densityA']
g = physical_parameters['gravity']
=======
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
g = pparams['gravity']
>>>>>>> TwoPhaseFlow

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
<<<<<<< HEAD
PRESSURE_model=3
V_model=1
PINC_model=2
=======
PRESSURE_model = mparams.pressure['index']
V_model = mparams.rans3p['index']
PINC_model = mparams.pressureIncrement['index']
>>>>>>> TwoPhaseFlow

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
<<<<<<< HEAD
dirichletConditions = {0: boundaryConditions['pressure_DBC']} 
advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC']}
=======
if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['pressure_DBC']}
    advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC']}
else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython()}
    advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython()}
>>>>>>> TwoPhaseFlow
