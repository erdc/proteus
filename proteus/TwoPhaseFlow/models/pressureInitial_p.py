from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import PresInit

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
<<<<<<< HEAD
=======

>>>>>>> TwoPhaseFlow
myTpFlowProblem = ct.myTpFlowProblem 
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# DOMAIN #
domain = myTpFlowProblem.domain

<<<<<<< HEAD
# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
PINIT_model=4
V_model=1
PRESSURE_model=3
=======
params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
PINIT_model = mparams.pressureInitial['index']
V_model = mparams.rans3p['index']
PRESSURE_model = mparams.pressure['index']
>>>>>>> TwoPhaseFlow

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
<<<<<<< HEAD
dirichletConditions = {0: boundaryConditions['pressure_DBC']}
advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC']}
diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['pressure_increment_DFBC']}}
=======
if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['pressure_DBC']}
    advectiveFluxBoundaryConditions = {0: boundaryConditions['pressure_AFBC']}
    diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['pressure_increment_DFBC']}}
else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].pInit_dirichlet.init_cython()}
    advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].pInit_advective.init_cython()}
    diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].pInit_diffusive.init_cython()}}
>>>>>>> TwoPhaseFlow
