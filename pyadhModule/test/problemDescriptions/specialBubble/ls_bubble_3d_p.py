from pyadh import *
from pyadh.default_p import *
from bubble3d import *

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

##\ingroup test
#\file ls_bubble_2d_p.py
#
# \todo finish ls_bubble_2d_p.py

coefficients = NCLevelSetCoefficients(V_model=1,RD_model=4,ME_model=2,VOF_model=3)

analyticalSolutions = None

#bubble rise
initialConditions  = {0:Bubble_phi(center=[bubbleCenter_x,bubbleCenter_y,bubbleCenter_z],radius=bubbleRadius)}

fluxBoundaryConditions = {0:'outFlow'}

dirichletConditions = {0:getDBC_phi}

diffusiveFluxBoundaryConditions = {0:{}}
