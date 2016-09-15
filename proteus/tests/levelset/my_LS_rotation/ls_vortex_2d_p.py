from proteus import *
from proteus.default_p import *
from math import *
from vortex2D import *
from proteus.mprans import NCLS
from proteus.ctransportCoefficients import smoothedHeaviside
#import Profiling

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

class OscillatingVortex2D(NCLS.Coefficients):
    #cek changed to put sphere inside arbitrary box with dimensions in L
    def __init__(self,L):
        self.radius = 0.25
        self.xc=0.0
        self.yc=0.5
    def uOfXT(self,x,t):
        return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
class OscillatingVortex2Dcylinder(NCLS.Coefficients):
    def uOfXT(self,x,t):
        L=1.0
        radius = 0.15*L
        xc=0.5*L
        yc=0.75*L
        init_condition = radius - math.sqrt((x[0]-xc)**2 + (x[1]-yc)**2)
        return init_condition
        #return smoothedHeaviside(epsFactHeaviside*he,init_condition)

analyticalSolution = {0:OscillatingVortex2Dcylinder()}

class UnitSquareVortex(TransportCoefficients.TC_base):
    from proteus.ctransportCoefficients import unitSquareVortexEvaluate
    from proteus.ctransportCoefficients import unitSquareVortexLevelSetEvaluate
    def __init__(self,useHJ=False,epsFact=1.5,checkMass=False,
                 RD_model=None,
                 useMetrics=0.0,sc_uref=1.0,sc_beta=1.0, 
                 cE=cE,cMax=cMax,ENTROPY_VISCOSITY=ENTROPY_VISCOSITY,SUPG=SUPG):
        self.waterline_interval=-1
        self.epsFact=epsFact
        self.useHJ = useHJ
        self.RD_modelIndex=RD_model
 	self.sc_uref=sc_uref
	self.sc_beta=sc_beta
	self.useMetrics=useMetrics
        mass={0:{0:'linear'}}
        advection={0:{0:'linear'}}
        diffusion={}
        potential={}
        reaction={}
        if self.useHJ:
            hamiltonian={0:{0:'linear'}}
        else:
            hamiltonian={}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
        self.checkMass=checkMass
        self.useMetrics = 0.0
	self.sc_uref=1.0
	self.sc_beta=1.0
        self.cE=cE
        self.cMax=cMax
        self.ENTROPY_VISCOSITY=ENTROPY_VISCOSITY
        self.SUPG=SUPG
    def attachModels(self,modelList):
        self.model = modelList[0]
	self.u_old_dof = numpy.copy(self.model.u[0].dof)
        self.q_v = numpy.zeros(self.model.q[('dH',0,0)].shape,'d')
        self.ebqe_v = numpy.zeros(self.model.ebqe[('dH',0,0)].shape,'d')
        self.unitSquareVortexLevelSetEvaluate(self.model.timeIntegration.tLast,
                                              self.model.q['x'],
                                              self.model.q[('u',0)],self.model.q[('grad(u)',0)],
                                              self.model.q[('m',0)],self.model.q[('dm',0,0)],
                                              self.model.q[('dH',0,0)],self.model.q[('dH',0,0)],
                                              self.model.q[('H',0)],self.q_v)
        self.model.q[('velocity',0)]=self.q_v
        self.model.ebqe[('velocity',0)]=self.ebqe_v
        if self.RD_modelIndex != None:
            #print self.RD_modelIndex,len(modelList)
            self.rdModel = modelList[self.RD_modelIndex]
        else:
            self.rdModel = self.model
    def preStep(self,t,firstStep=False):
        pi = math.pi
        x = self.model.q['x'][...,0]
        y = self.model.q['x'][...,1]
        self.q_v[...,0]  = -2.0*pi*(y-0.5)
        self.q_v[...,1]  =  2.0*pi*(x-0.5)
        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
       	self.u_old_dof = numpy.copy(self.model.u[0].dof)
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        pass

coefficients = UnitSquareVortex(useHJ=True,epsFact=epsFactHeaviside,checkMass=checkMass,useMetrics=useMetrics)

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    #if (x[1] == 0.0):
    #    return lambda x,t: 0.0
    #if (x[0] == 0.0 or
    #    x[0] == 1.0 or
    #    x[1] == 0.0 or
    #    x[1] == 1.0):
    #    return lambda x,t: 0.0
def zeroInflow(x):
    return lambda x,t: 0.0
    # if (x[0] == 0.0 and x[1] <= 0.5):
    #     return lambda x,t: 0.0
    # if (x[0] == 1.0 and x[1] >= 0.5):
    #     return lambda x,t: 0.0
    # if (x[1] == 0.0 and x[0] >= 0.5):
    #     return lambda x,t: 0.0
    # if (x[1] == 1.0 and x[0] <= 0.5):
    #     return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x):
    return lambda x,t: 0.0
advectiveFluxBoundaryConditions =  {}
#advectiveFluxBoundaryConditions =  {0:zeroadv}


diffusiveFluxBoundaryConditions = {0:{}}

## @}
