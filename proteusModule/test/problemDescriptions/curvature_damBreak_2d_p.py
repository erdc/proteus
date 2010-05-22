from pyadh import *
from pyadh.default_p import *
from damBreak import *
from math import *

coefficients = LevelSetCurvatureCoefficients(epsFact=epsFact_curvature,
                                             LSModel_index=2,
                                             nd=nd)
movingDomain=True
fluxBoundaryConditions = {0:'outFlow'}

def getDBC_kappa(x,flag):
    pass

dirichletConditions = {0:getDBC_kappa}

def getAFBC_kappa(x,flag):
     pass

advectiveFluxBoundaryConditions =  {0:getAFBC_kappa}

diffusiveFluxBoundaryConditions = {0:{}}
