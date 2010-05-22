from pyadh import *
from pyadh.default_p import *
import sys
sys.path.append("..")
import LN_Example_Domains
from pyadh import LADRellam
"""
linear advection dispersion reaction example in
flow field from Mario Putti's Larson and Niklasson example 1

Domain is heterogeneous with 2 low conductivity blocks 
"""
name = "LN_Ex1_adr"
##\page Tests Test Problems 
# \ref LN_Example1_adr_p.py "Linear ADR in Heterogeneous domain"
#

##\ingroup test
#\file LN_Example1_adr_p.py
#
#\brief linear advection-dispersion in single phase flow in block heterogeneous domain
#constant head on left and right
from LN_Example1 import *
if useELLAM:
    LevelModelType = LADRellam.OneLevelLADR

coefficients = GroundwaterTransportCoefficientsELLAM(omega=omega,
                                                     alpha_L=alpha_L,
                                                     alpha_T=alpha_T,
                                                     d = d,
                                                     nd = nd,
                                                     velocitySpaceFlag = velocitySpaceFlag,
                                                     meModelId  =1,
                                                     flowModelId=0)


def inflowFlux(x,t):
    if x[1] <= yInflowStart:
        return 0.0
    if (smoothInflow and t <= smoothDt):
        return inflowVal*(t-0.0)/smoothDt
    return inflowVal

    

def getDBC(x,tag):
    pass
dirichletConditions = {0:getDBC}

analyticalSolution = analyticalSolution_trans

initialConditions  = initialConditions_trans

fluxBoundaryConditions = {0:'outFlow'}#'noFlow'}

def getAFBC(x,tag):
    if x[0] == 0.0:
        return lambda x,t: inflowFlux(x,t)

advectiveFluxBoundaryConditions =  {0:getAFBC}


def getDFBC(x,tag):
   if x[1] in [0.0,1.0]: 
       return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}


