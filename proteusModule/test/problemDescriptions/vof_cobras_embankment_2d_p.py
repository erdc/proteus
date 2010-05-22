from pyadh import *
from pyadh.default_p import *
from embankmentRunup import *

coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2)


dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Flat_H()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
