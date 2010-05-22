from pyadh import *
from pyadh.default_p import *
from susan import *
"""
The non-conservative level set description of the free surface in two-phase flow around the bow of the SUSAN MAERSK in 2D
"""

##

##\ingroup test
#\file ls_susan_2d_p.py
#
#\brief The non-conservative level set description of the free surface in two-phase flow around the bow of the SUSAN MAERSK in 2D

analyticalSolutions = None

coefficients = NCLevelSetCoefficients()

def getDBC(x):
    if x[0] >= (2.0e2-1.0e-8):
        return lambda x,t: x[1] - waterLevel

dirichletConditions = {0:getDBC}

class PerturbedSurface_phi:
    def __init__(self,waterLevel,slope):
        self.waterLevel=waterLevel
        self.slope=slope
    def uOfXT(self,x,t):
        z = self.waterLevel + (x[0] - 0.5*L[0])*self.slope
        return x[1]-z

initialConditions  = {0:PerturbedSurface_phi(waterLevel,slope)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
