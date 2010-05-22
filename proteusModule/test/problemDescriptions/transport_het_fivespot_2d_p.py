from pyadh import *
from pyadh.default_p import *
from gw_het_fivespot import *

"""
advective-dispersive transport in heterogeneous domain
"""

##\page Tests Test Problems 
# \ref transport_het_fivespot_2d_p.py "linear transport equation in heterogeneous porous medium"
#

##\ingroup test
#\file transport_het_fivespot_2d_p.py
#
#\brief linear transport equation in heterogeneous porous medium

initialConditions = {0:constantIC(backgroundConc)}

analyticalSolution = transportAnalyticalSolution

dirichletConditions = {0:inflowConcentration}

coefficients = GroundwaterTransportCoefficients(nctrans,
                                                omega=omega,
                                                alpha_L=alpha_L,
                                                alpha_T=alpha_T,
                                                d = d_mol)
   


fluxBoundaryConditions = {0:'outFlow'}
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}


