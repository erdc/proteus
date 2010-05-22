from pyadh import *
from pyadh.default_p import *
from wavetank import *

if applyCorrection and applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1)
elif applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=2,ME_model=1)
else:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=None,ME_model=1)
    
dirichletConditions = {0:getDBC_phi}

initialConditions  = {0:Flat_phi()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

