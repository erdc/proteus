from pyadh import *
from pyadh.default_p import *
from simpleweir import *
"""
The non-conservative level-set representation of the free surface in two-phase flow over a weir.
"""

##

##\ingroup test
#\file ls_so_simpleweir_2d_p.py
#\brief The non-conservative level-set representation of the free surface in two-phase flow over a weir.
#\todo finish ls_so_simpleweir_2d_p.py doc
#\todo figure out if doxygen can extract some of the information so we don't have to repeat so much in the comments
#\todo put commonly used equation sets in doc.h so we can link to them in the p-file doc
analyticalSolutions = None

coefficients = NCLevelSetCoefficients()

def getDBC(x):
    pass
#     if x[0] == 0.0 or x[0] >= flumeEnd - 1.0e-8:
#         return lambda x,t: x[1]-waterLevel

dirichletConditions = {0:getDBC}

class StraightIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return x[1] - waterLevel

initialConditions  = {0:StraightIC()}

#shouldn't matter for noncnsrv formulation
fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}


