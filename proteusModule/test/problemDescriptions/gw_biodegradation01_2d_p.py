from pyadh import *
from pyadh.default_p import *
from gw_biodegradation01 import *

"""
Heterogenous Poisson's equation in 2D
"""

##\page Tests Test Problems 
# \ref"Heterogeneous Poisson's equation"
#

##\ingroup test
#\file 
#
#\brief Heterogeneous Poisson's equation in 2D

initialConditions = None

Ident = numpy.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0

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



##################################################


aOfX = {0:hydraulicConductivity}; fOfX = {0:fluidSourceTerms}

analyticalSolution = flowAnalyticalSolution
dirichletConditions = {0:flowHeadBCs}

analyticalSolutionVelocity = velocityAnalyticalSolution

coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,ncflow)
   

fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {0:dummyBCs}

diffusiveFluxBoundaryConditions = {0:{0:noFlowBCs}}


