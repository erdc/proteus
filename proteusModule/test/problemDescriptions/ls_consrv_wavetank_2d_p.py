from pyadh import *
from pyadh.default_p import *
from wavetank import *

coefficients = LevelSetConservation(applyCorrection=True,
                                    epsFactHeaviside=epsFact_massCorrection_heaviside,
                                    epsFactDirac=epsFact_massCorrection_dirac,
                                    epsFactDiffusion=epsFact_massCorrection_diffusion,
                                    LSModel_index=1,V_model=0,me_model=4,VOFModel_index=2,nd=nd)

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
