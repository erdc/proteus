from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import MoveMeshMonitor
import numpy as np
from proteus import Context

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
mparams = myTpFlowProblem.Models # model parameters
pparams = myTpFlowProblem.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
nmediatypes = len(domain.regionflags)  # (!) should be region flags
smtypes     = np.zeros((nmediatypes+1, 2), 'd')
smflags     = np.zeros((nmediatypes+1,), 'i')
smtypes[:, 0] = 1.0    # e
smtypes[:, 1] = 0.3    # nu

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = mparams.movemeshmonitor['index']
assert ME_model is not None, 'vof model index was not set!'
if mparams.movemeshmonitor['useLS'] is True:
    LS_MODEL = mparams.ncls['index']
else:
    LS_MODEL = None
initialConditions = None

analyticalSolution = {}

LevelModelType = MoveMesh.LevelModel
coefficients = MoveMesh.Coefficients(ndnd,
                                     V_model=V_model,
                                     modelType_block=smFlags,
                                     modelParams_block=smTypes,
                                     meIndex=ME_model,
                                     LS_MODEL=LS_MODEL)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = None

if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['hx'],
                           1: boundaryConditions['hy']}
    stressFluxBoundaryConditions = {0: boundaryConditions['u_stress'],
                                    1: boundaryConditions['v_stress']}
    if nd == 3:
        dirichletConditions[2] = boundaryConditions['hz']
        stressFluxBoundaryConditions[2] = boundaryConditions['w_stress']

else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].hx_dirichlet.init_cython(),
                           1: lambda x, flag: domain.bc[flag].hy_dirichlet.init_cython()}
    stressFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].u_stress,
                                    1: lambda x, flag: domain.bc[flag].v_stress}

fluxBoundaryConditions = {0: 'noFlow',
                          1: 'noFlow'}

advectiveFluxBoundaryConditions = {}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {}}

if nd == 3:
    fluxBoundaryConditions[2] = 'noFlow'
    diffusiveFluxBoundaryConditions[2] = {}

