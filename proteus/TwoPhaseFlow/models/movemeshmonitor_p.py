from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import MoveMeshMonitor
import numpy as np
from proteus import Context

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions
params = ct.params.movemeshmonitor

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
func = params['func']
he_max = params['he_max']
he_max = params['he_min']
ntimes_solved = params['ntimes_solved']
fixedNodeMaterialTypes = params['fixedNodeMaterialTypes']
fixedElementMaterialTypes = params['fixedElementMaterialTypes']
noNodeVelocityMaterialTypes = params['noNodeVelocityMaterialTypes']
nSmoothOut = params['nSmoothOut']
epsFact = params['epsFact']
epsTimeStep = params['epsTimeStep']
grading = params['grading']
grading_type = params['grading_type']
resetNodeVelocityArray = params['resetNodeVelocityArray']
useLS = params['useLS']
do_firstStep = params['do_firstStep']
func = params['func']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = params['index']
assert ME_model != None, 'vof model index was not set!'
if params['useLS'] is True:
    LS_MODEL = params.ncls['index']
else:
    LS_MODEL = False
initialConditions = None

analyticalSolution = {}

coefficients = MoveMeshMonitor.Coefficients(func=func,
                                            nd=nd,
                                            he_max=he_max,
                                            he_min=he_min,
                                            ntimes_solved=movemesh_ntimes,
                                            LS_MODEL=LS_MODEL,
                                            ME_MODEL=ME_MODEL,
                                            fixedNodeMaterialTypes=fixedNodeMaterialTypes,
                                            fixedElementMaterialTypes=fixedElementMaterialTypes,
                                            noNodeVelocityNodeMaterialTypes=noNodeVelocityNodeMaterialTypes,
                                            nSmoothOut=nSmoothOut,
                                            epsTimeStep=epsTimeStep,
                                            epsFact_density=epsFact,
                                            grading=grading,
                                            grading_type=grading_type,
                                            scale_with_nd=scale_with_nd,
                                            do_firstStep=do_firstStep)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = None

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
def getDBC(x, flag):
    if flag != 0:
        return None
dirichletConditions = {0: getDBC}

def getDFBC(x, flag):
    return lambda x, t: 0.
diffusiveFluxBoundaryConditions = {0: {0: getDFBC}}
