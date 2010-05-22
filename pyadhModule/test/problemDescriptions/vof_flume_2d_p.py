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

coefficients = VOFCoefficients(LS_model=2,V_model=1,RD_model=4,ME_model=3)

analyticalSolutions = None
dirichletConditions = {0:getDBC_vof}
initialConditions  = {0:Flat_H()}
fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
