from pyadh import *
from pyadh.default_p import *
from wavetank3d import *
from math import *

coefficients = LevelSetCurvatureCoefficients(epsFact=epsFact_curvature,LSModel_index=2,nd=nd)

fluxBoundaryConditions = {0:'outFlow'}

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

def getAFBC(x,flag):
     pass

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
