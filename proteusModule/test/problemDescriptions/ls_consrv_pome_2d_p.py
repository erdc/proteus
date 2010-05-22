from pyadh import *
from pyadh.default_p import *
from pome import *

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

##\ingroup test
#\file ls_bubble_2d_p.py
#
# \todo finish ls_bubble_2d_p.py

coefficients = LevelSetConservation(LSModel_index=3,V_model=2,me_model=5,VOFModel_index=4)

waterLevel = 0.9*L[1]

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
#     if x[0] == 0.0:
#         if x[1] == 0.0:
#             return lambda x,t: 0.0


dirichletConditions = {0:getDBC_cnsrv}
#bubble rise
initialConditions  = {0:zero_phi()}

fluxBoundaryConditions = {0:'outFlow'}


def getAFBC_cnsrv(x):
    pass
#     #should really be plain outflow but want to ensure it's zero for bubble
#     if (x[0] == 0.0 or
#         x[0] == L[0] or
#         x[1] == 0.0 or
#         x[1] == L[1]):
#         return lambda x,t: 0.0
def getDFBC_cnsrv(x):
    #should really be plain outflow but want to ensure it's zero for bubble
    if (x[0] == 0.0 or
        x[0] == L[0] or
        x[1] == 0.0 or
        x[1] == L[1]):
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_cnsrv}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_cnsrv}}
