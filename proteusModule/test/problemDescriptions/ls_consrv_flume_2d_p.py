from pyadh import *
from pyadh.default_p import *
from flume import *

"""
The non-conservative level set description of a flume in a two-phase flow
"""

##\ingroup test
#\file ls_flume_2d_p.py
#
# \todo finish ls_flume_2d_p.py

coefficients = LevelSetConservation(applyCorrection=applyCorrection,
                                    epsFactHeaviside=epsFact_massCorrection_heaviside,
                                    epsFactDirac=epsFact_massCorrection_dirac,
                                    epsFactDiffusion=epsFact_massCorrection_diffusion,
                                    LSModel_index=2,V_model=1,me_model=5,VOFModel_index=3,nd=nd)

class zero_phi:
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return self.uOfX(X)

analyticalSolutions = None

def getDBC_cnsrv(x,flag):
    pass


dirichletConditions = {0:getDBC_cnsrv}
initialConditions  = {0:zero_phi()}

fluxBoundaryConditions = {0:'outFlow'}


def getAFBC_cnsrv(x,flag):
    pass

def getDFBC_cnsrv(x,flag):
    pass 

advectiveFluxBoundaryConditions =  {0:getAFBC_cnsrv}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_cnsrv}}
