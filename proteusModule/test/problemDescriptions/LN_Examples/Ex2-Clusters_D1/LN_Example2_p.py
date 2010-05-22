from pyadh import *
from pyadh.default_p import *


"""
Heterogeneous (2 block) Poisson equation

Mario Putti's Larson and Niklasson example 2
"""
name = "LN_Example2_p1_pwl"
##\page Tests Test Problems 
# \ref LN_Example1_p.py "Heterogeneous Poisson's equation"
#

##\ingroup test
#\file LN_Example1_p.py
#
#\brief single phase flow in block heterogeneous domain
#constant head on left and right
from LN_Example2 import *
eps = 1.0e-6
def headBCs(x,tag):
    if x[0] <= 1.0e-6:
        return lambda x,t: head_left
    if x[0] >= domain.L[0]-eps:
        return lambda x,t: head_right
def noFlowBCs(x,tag):
    if x[1] <= 1.0e-6 or x[1] >= domain.L[1]-eps:
        return lambda x,t: 0.0
    
class velEx:
    def __init__(self,duex,aex):
        self.duex = duex
        self.aex = aex
    def uOfX(self,X):
        du = self.duex.duOfX(X)
        A  = numpy.reshape(self.aex(X),(2,2))
        return -numpy.dot(A,du)
    def uOfXT(self,X,T):
        return self.uOfX(X)


initialConditions = initialConditions_flow


analyticalSolution = None
dirichletConditions = {0:headBCs}

analyticalSolutionVelocity = None

coefficients = SubsurfaceTransportCoefficients.SinglePhaseDarcyCoefficients(hydraulicConductivities,sources,nc=1,nd=nd,
                                                                            materialValuesLocallyConstant=True,
                                                                            timeVaryingCoefficients=False)
   

fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {0:noFlowBCs}#{0:dummyBCs}

diffusiveFluxBoundaryConditions = {0:{0:noFlowBCs}}


