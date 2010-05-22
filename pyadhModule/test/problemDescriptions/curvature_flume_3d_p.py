from pyadh import *
from pyadh.default_p import *
from flume3d import *
from math import *

coefficients = LevelSetCurvatureCoefficients(epsFact=epsFact_curvature,LSModel_index=2)

fluxBoundaryConditions = {0:'outFlow'}
def getDBC(x):
    pass

dirichletConditions = {0:getDBC}

def getAFBC(x):
     pass
advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
