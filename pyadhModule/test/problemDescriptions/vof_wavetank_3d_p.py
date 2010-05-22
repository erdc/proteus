from pyadh import *
from pyadh.default_p import *
from wavetank3d import *

coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2)

analyticalSolutions = None

dirichletConditions = {0:getDBC_vof_wavetank}

initialConditions  = {0:Flat_H()}

fluxBoundaryConditions = {0:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof_wavetank}

diffusiveFluxBoundaryConditions = {0:{}}
