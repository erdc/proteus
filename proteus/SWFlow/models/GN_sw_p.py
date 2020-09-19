from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import GN_SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
mySWFlowProblem = ct.mySWFlowProblem
genMesh=mySWFlowProblem.genMesh
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
LevelModelType = GN_SW2DCV.LevelModel
coefficients = GN_SW2DCV.Coefficients(g=g,
                                      bathymetry={
                                          0: bathymetry} if bathymetry is not None else None,
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
                     4: initialConditions['h_times_w'],
                     5: initialConditions['h_times_beta']}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: boundaryConditions['water_height'],
                       1: boundaryConditions['x_mom'],
                       2: boundaryConditions['y_mom'],
                       3: boundaryConditions['h_times_eta'],
                       4: boundaryConditions['h_times_w'],
                       5: boundaryConditions['h_times_beta']}
fluxBoundaryConditions = {0: 'outFlow',
                          1: 'outFlow',
                          2: 'outFlow',
                          3: 'outFlow',
                          4: 'outFlow',
                          5: 'outFlow'}
advectiveFluxBoundaryConditions = {0: lambda x, flag: None,
                                   1: lambda x, flag: None,
                                   2: lambda x, flag: None,
                                   3: lambda x, flag: None,
                                   4: lambda x, flag: None,
                                   5: lambda x, flag: None}
diffusiveFluxBoundaryConditions = {0: {},
                                   1: {1: lambda x, flag: None},
                                   2: {2: lambda x, flag: None},
                                   3: {3: lambda x, flag: None},
                                   4: {4: lambda x, flag: None},
                                   5: {5: lambda x, flag: None}}
# **************************************** #
# ********** ANALYTICAL SOLUTION ********* #
# **************************************** #
if (mySWFlowProblem.analyticalSolution is not None):
    analyticalSolution = {0: analyticalSolution['h_exact'],
                          1: analyticalSolution['hu_exact'],
                          2: analyticalSolution['hv_exact'],
                          3: analyticalSolution['heta_exact'],
                          4: analyticalSolution['hw_exact'],
                          5: analyticalSolution['hbeta_exact']}
