from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import MCorr
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
initialConditions = myTpFlowProblem.initialConditions
boundaryConditions = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain

# DOMAIN #
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
myparams = params.Models.mcorr
pparams = params.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = mparams.mcorr.index
assert ME_model is not None, 'vof model index was not set!'
LS_model = mparams.ncls.index
VOF_model = mparams.vof.index
if mparams.rans2p.index is not None:
    V_model = mparams.rans2p.index
elif mparams.rans3p.index is not None:
    V_model = mparams.rans3p.index
else:
    assert mparams.rans2p.index is not None or params.rans3p.index is not None, 'RANS2P or RANS3P must be used with VOF'

LevelModelType = MCorr.LevelModel

coefficients = MCorr.Coefficients(LSModel_index=LS_model,
                                  V_model=V_model,
                                  me_model=ME_model,
                                  VOFModel_index=VOF_model,
                                  applyCorrection=myparams.applyCorrection,
                                  nd=nd,
                                  checkMass=myparams.checkMass,
                                  useMetrics=myparams.useMetrics,
                                  epsFactHeaviside=myparams.epsFactHeaviside,
                                  epsFactDirac=myparams.epsFactDirac,
                                  epsFactDiffusion=myparams.epsFactDiffusion)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
class zero_phi:
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return 0.0
initialConditions  = {0:zero_phi()}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
# N/A
