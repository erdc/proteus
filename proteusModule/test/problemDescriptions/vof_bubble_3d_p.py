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

coefficients = VOFCoefficients(LS_model=2,V_model=1,RD_model=4,ME_model=3)


analyticalSolutions = None


dirichletConditions = {0:getDBC_vof}

#bubble rise
initialConditions  = {0:Bubble_H(center=[bubbleCenter_x,bubbleCenter_y,bubbleCenter_z],radius=bubbleRadius)}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
