from pyadh import *
from pyadh.default_p import *
"""
The non-conservative level set description of the free surface in two-phase flow over a weir.
"""

##

##\ingroup test
#\file ls_so_weir_2d_p.py
#
#\brief The non-conservative level set description of the free surface in two-phase flow over a weir.

nd = 2

analyticalSolutions = None

polyfile = "simpleWeir"

coefficients = NCLevelSetCoefficients()

def getDBC(x):
    if x[0] == 1.0:
        return lambda x,t:  x[1] - 0.5

dirichletConditions = {0:getDBC}

class StraightIC:
    def __init__(self,dbc):
        self.dbc=dbc
    def uOfXT(self,x,t):
        return self.dbc([1.0,0.0,0.0])(x,t)

initialConditions  = {0:StraightIC(getDBC)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 100
