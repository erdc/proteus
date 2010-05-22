from pyadh import *
from pyadh.default_p import *
from math import *
from wavetank import *

"""
The redistancing model used by the wavetank test problem.
"""
##

##\ingroup test
#\file redist_wavetank_2d_p.py
#\brief The redistancing model used by the wavetank test problem.

if applyCorrection and applyRedistancing:
    coefficients = RedistanceLevelSet(applyRedistancing=True,
                                      epsFact=epsFact_redistance,
                                      nModelId=1,
                                      rdModelId=3)
else:
    coefficients = RedistanceLevelSet(applyRedistancing=True,
                                      epsFact=epsFact_redistance,
                                      nModelId=1,
                                      rdModelId=2)

#now define the Dirichlet boundary conditions

dirichletConditions = {0:getDBC_phi}

weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
