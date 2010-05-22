from pyadh import *
from pyadh.default_p import *
from porousObstacle import *

"""
The non-conservative level set description of a porousObstacle in a two-phase flow
"""

##\ingroup test
#\file ls_porousObstacle_2d_p.py
#
# \todo finish ls_porousObstacle_2d_p.py
name = "LevelSet"

coefficients = NCLevelSetCoefficients(V_model=1,RD_model=4,ME_model=2,VOF_model=3)

analyticalSolutions = None

dirichletConditions = {0:getDBC_phi}
initialConditions  = {0:Flat_phi()}
fluxBoundaryConditions = {0:'outFlow'}
