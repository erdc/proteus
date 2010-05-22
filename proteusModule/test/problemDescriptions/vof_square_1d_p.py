from pyadh import *
from pyadh.default_p import *
from math import *
from square import *

class VOFSquareWave1D:
    def __init__(self):
        self.radius = 0.25
        self.velocity=1.0
    def uOfXT(self,x,t):
        center = (0.5 + t)%1.0
        phi = self.radius - abs(x[0] - center)
        if phi > 0:
            return 1.0
        elif phi < 0:
            return 0.0
        else:
            return 0.5

analyticalSolution = {0:VOFSquareWave1D()}

class SquareWave1D(TransportCoefficients.VOFCoefficients):
    def __init__(self,LS_model=-1,RD_model=-1,ME_model=1,EikonalSolverFlag=0,checkMass=True,epsFact=0.0):
        TransportCoefficients.VOFCoefficients.__init__(self,
                                                       LS_model=LS_model,
                                                       V_model=-1,
                                                       RD_model=RD_model,
                                                       ME_model=ME_model,
                                                       EikonalSolverFlag=EikonalSolverFlag,
                                                       checkMass=checkMass,
                                                       epsFact=epsFact)
    def initializeElementQuadrature(self,t,cq):
        self.q_v = numpy.ones(cq[('f',0)].shape,'d')
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.ebq_v = numpy.ones(cebq[('f',0)].shape,'d')
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        self.ebqe_v = numpy.ones(cebqe[('f',0)].shape,'d')

if applyRedistancing:
    coefficients = SquareWave1D(LS_model=0,RD_model=1,ME_model=2)
else:
    coefficients = SquareWave1D(LS_model=0,ME_model=1)
    
coefficients.variableNames=['vof']

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    
dirichletConditions = {0:getDBC}

eps=1.0e-8
def getPDBC(x,flag):
    if x[0] <= eps or x[0] >= (1.0-eps):
        return numpy.array([0.0,0.0,0.0])

#periodic boundary conditions

periodicDirichletConditions = {0:getPDBC}

initialConditions  = {0:analyticalSolution[0]}


def getAFBC(x,flag):
    pass
#     if x[0] == 0.0:
#         return lambda x,t: 0.0
#     if x[0] == 1.0:
#         return lambda x,t: 0.0
    
fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}

## @}
