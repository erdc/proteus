from pyadh import *
from pyadh.default_p import *
from math import *
from square import *
from pyadh import NCLS
if useNCLS:
    LevelModelType = NCLS.OneLevelNCLS
class SquareWave1D:
    def __init__(self):
        self.radius = 0.25
        self.velocity=1.0
    def uOfXT(self,x,t):
        center = (0.5 + t)%1.0
#        return math.cos(math.pi*(x[0] - center))**2
        return self.radius - abs(x[0] - center)

analyticalSolution = {0:SquareWave1D()}

class SquareWave1D(TransportCoefficients.NCLevelSetCoefficients):
    def __init__(self,
                 RD_model=-1,
                 ME_model=1,
                 EikonalSolverFlag=0,
                 checkMass=True):
        TransportCoefficients.NCLevelSetCoefficients.__init__(self,
                                                              V_model=-1,
                                                              RD_model=RD_model,
                                                              ME_model=ME_model,
                                                              EikonalSolverFlag=EikonalSolverFlag,
                                                              checkMass=checkMass)
    def initializeElementQuadrature(self,t,cq):
        self.q_v = numpy.ones(cq[('grad(u)',0)].shape,'d')
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.ebq_v = numpy.ones(cebq[('grad(u)',0)].shape,'d')
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        self.ebqe_v = numpy.ones(cebqe[('grad(u)',0)].shape,'d')

if applyRedistancing:
    coefficients = SquareWave1D(RD_model=1,
                                ME_model=0,
                                EikonalSolverFlag=fmmFlag,
                                checkMass=checkMass)
else:
    coefficients = SquareWave1D(RD_model=-1,
                                ME_model=0,
                                EikonalSolverFlag=fmmFlag,
                                checkMass=checkMass)
    

coefficients.variableNames=['phi']

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    
dirichletConditions = {0:getDBC}

def getPDBC(x,flag):
    if x[0] == 0.0 or x[0] == 1.0:
        return numpy.array([0.0,0.0,0.0])

#periodic boundary conditions

periodicDirichletConditions = {0:getPDBC}

initialConditions  = {0:analyticalSolution[0]}

advectiveFluxBoundaryConditions =  {}#can't have for dg {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}

## @}
