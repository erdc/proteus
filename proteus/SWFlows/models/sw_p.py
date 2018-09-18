from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import Context

# READ FROM CONTEXT #
ct = Context.get()
initialConditions = ct.initialConditions
boundaryConditions = ct.boundaryConditions
bathymetry_function = ct.bathymetry_function

# PHYSICAL PARAMETERS #
g=ct.physical_parameters['g']
LINEAR_FRICTION=ct.physical_parameters['LINEAR_FRICTION']
mannings=ct.physical_parameters['mannings']

# NUMERICAL PARAMETERS #
cE=ct.numerical_parameters['cE']
LUMPED_MASS_MATRIX=ct.numerical_parameters['LUMPED_MASS_MATRIX']

# DOMAIN #
nd = 2
domain = ct.domain

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: initialConditions['water_height'],
                     1: initialConditions['x-mom'],
                     2: initialConditions['y-mom']}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
dirichletConditions = {0: boundaryConditions['water_height'],
                       1: boundaryConditions['x-mom'],
                       2: boundaryConditions['y-mom']}

fluxBoundaryConditions = {0: 'outFlow',
                          1: 'outFlow',
                          2: 'outFlow'}
advectiveFluxBoundaryConditions =  {0: lambda x,flag: lambda x,t: 0.0,
                                    1: lambda x,flag: lambda x,t: 0.0,
                                    2: lambda x,flag: lambda x,t: 0.0}
diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1: lambda x,flag: lambda x,t: 0.0},
                                   2:{2: lambda x,flag: lambda x,t: 0.0}}

#########################################
##### CREATE MODEL AND COEFFICIENTS #####
#########################################
bathymetry={0: bathymetry_function}
LevelModelType = SW2DCV.LevelModel
coefficients = SW2DCV.Coefficients(g=g,
                                   bathymetry={0: bathymetry_function},
                                   cE=cE,
                                   LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,
                                   LINEAR_FRICTION=LINEAR_FRICTION,
                                   mannings=mannings)
