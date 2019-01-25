from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import DSW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
mySWFlowProblem = ct.mySWFlowProblem
physical_parameters = mySWFlowProblem.physical_parameters
numerical_parameters = mySWFlowProblem.swe_parameters
initialConditions = mySWFlowProblem.initialConditions
boundaryConditions = mySWFlowProblem.boundaryConditions
bathymetry = mySWFlowProblem.bathymetry
reflecting_BCs = mySWFlowProblem.reflectingBCs
analyticalSolution = mySWFlowProblem.analyticalSolution

# DOMAIN #
nd = 2
domain = mySWFlowProblem.domain
if domain is None:
    meshfile = mySWFlowProblem.AdH_file

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
# PHYSICAL PARAMETERS #
g = physical_parameters['gravity']
LINEAR_FRICTION = physical_parameters['LINEAR_FRICTION']
mannings = physical_parameters['mannings']

# NUMERICAL PARAMETERS #
cE = numerical_parameters['cE']
LUMPED_MASS_MATRIX = numerical_parameters['LUMPED_MASS_MATRIX']

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = DSW2DCV.LevelModel
coefficients = DSW2DCV.Coefficients(g=g,
                                   bathymetry={0:bathymetry} if bathymetry is not None else None,
                                   cE=cE,
                                   LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,
                                   LINEAR_FRICTION=LINEAR_FRICTION,
                                   mannings=mannings)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: initialConditions['water_height'],
                     1: initialConditions['x_mom'],
                     2: initialConditions['y_mom'],
                     3: initialConditions['h_times_eta'],
                     4: initialConditions['h_times_w']}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: boundaryConditions['water_height'],
                       1: boundaryConditions['x_mom'],
                       2: boundaryConditions['y_mom'],
                       3: boundaryConditions['h_times_eta'],
                       4: boundaryConditions['h_times_w']}
fluxBoundaryConditions = {0: 'outFlow',
                          1: 'outFlow',
                          2: 'outFlow',
                          3: 'outFlow',
                          4: 'outFlow'}
advectiveFluxBoundaryConditions =  {0: lambda x, flag: None,
                                    1: lambda x, flag: None,
                                    2: lambda x, flag: None,
                                    3: lambda x, flag: None,
                                    4: lambda x, flag: None}
diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1: lambda x,flag: None},
                                   2:{2: lambda x,flag: None},
                                   3:{3: lambda x,flag: None},
                                   4:{4: lambda x,flag: None}}
