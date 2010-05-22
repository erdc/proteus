from pyadh import *
from pyadh.default_p import *
from vortex import *



coefficients = LevelSetConservation(applyCorrection=applyCorrection,
#                                    epsFactHeaviside=0.0,#1.0e-2,
                                    epsFactHeaviside=1.0e-2,
                                    epsFactDirac=1.0e-2,
                                    epsFactDiffusion=1.0e2,
                                    LSModel_index=0,
                                    V_model=0,
                                    me_model=3,
                                    VOFModel_index=1)
# coefficients = LevelSetConservation(applyCorrection=applyCorrection,
#                                     epsFactHeaviside=0.0,
#                                     epsFactDirac=0.1,
#                                     epsFactDiffusion=1.0e2,
#                                     LSModel_index=0,
#                                     V_model=0,
#                                     me_model=3,
#                                     VOFModel_index=1)

class zero_phi:
    def __init__(self):
        pass
    def uOfX(self,X):
        return 0.0
    def uOfXT(self,X,t):
        return self.uOfX(X)

analyticalSolutions = None

def getDBC_cnsrv(x):
    pass

dirichletConditions = {0:getDBC_cnsrv}

initialConditions  = {0:zero_phi()}

fluxBoundaryConditions = {0:'noFlow'}


def getAFBC_cnsrv(x):
    return lambda x,t: 0.0
def getDFBC_cnsrv(x):
    return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_cnsrv}}
