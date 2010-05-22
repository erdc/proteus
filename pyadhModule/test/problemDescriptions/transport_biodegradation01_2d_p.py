from pyadh import *
from pyadh.default_p import *
from gw_biodegradation01 import *

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

initialConditions = {0:constantIC(backgroundConc[0]),1:constantIC(backgroundConc[1]),2:constantIC(backgroundConc[2])}

analyticalSolution = transportAnalyticalSolution

dirichletConditions = {0:inflowConcentration0,1:inflowConcentration1,2:inflowConcentration2}

coefficients = GroundwaterBiodegradation01Coefficients(omega=omega,
                                                       alpha_L=alpha_L,
                                                       alpha_T=alpha_T,
                                                       Kox_max=Kox_max,
                                                       Kox_C=Kox_C,
                                                       Kox_E=Kox_E,
                                                       Kox_X=Kox_X,
                                                       Yield=Yield,
                                                       k_d  =k_d)
   


fluxBoundaryConditions = {0:'outFlow',1:'outFlow',2:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}

#advectiveFluxBoundaryConditions =  {0:zeroInflow_tran}

#diffusiveFluxBoundaryConditions = {0:{0:zeroInflow_tran_diff}}

