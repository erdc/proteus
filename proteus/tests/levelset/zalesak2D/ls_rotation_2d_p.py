from proteus import *
from proteus.default_p import *
from math import *
from rotation2D import *
from proteus.mprans import NCLS

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

class Zalesak2D:
    def __init__(self,
                 L=[1.0,
                    1.0],
                 center=[0.5,
                         0.5],
                 radius=0.45,
                 zalesak=True):
        self.zalesak=zalesak
        self.radius = radius
        self.slotwidth = radius*3.0/9.0
        self.slotlength = radius
        self.xc = center
        self.xnw = [center[0] - 0.5*self.slotwidth,
                    center[1] - (radius - self.slotlength)]
        self.xne = [center[0] + 0.5*self.slotwidth,
                    center[1] - (radius - self.slotlength)]
        self.xsw = [center[0] - 0.5*self.slotwidth,
                    center[1] - (radius)]
        self.xse = [center[0] + 0.5*self.slotwidth,
                    center[1] - (radius)]
    def uOfXT(self,x,t):
        if not self.zalesak:
            return self.radius - math.sqrt((x[0]-self.xc[0])**2 + (x[1]-self.xc[1])**2)
        else:
            from math import sqrt
            dist = lambda u,v: sqrt( (u[0] - v[0])**2 + (u[1] - v[1])**2)
            phic = dist(self.xc,x) - self.radius
            phine = -dist(self.xne,x)
            phinw = -dist(self.xnw,x)
            phise = dist(self.xse,x)
            phisw = dist(self.xsw,x)
            phin = self.xnw[1] - x[1]
            phis = -(self.xsw[1] - x[1])
            phie = self.xne[0] - x[0]
            phiw = -(self.xnw[0] - x[0])
            if x[1] >= self.xnw[1]:
                if x[0] < self.xnw[0]:
                    phi = max(phic,phinw)
                else:
                    if x[0] < self.xne[0]:
                        phi = max(phic,phin)
                    else:
                        phi = max(phic,phine)
            elif x[1] >= self.xsw[1]:
                if x[0] < self.xnw[0]:
                    phi = max(phic,phiw)
                else:
                    if x[0] < self.xne[0]:
                        phi = min([phin,phie,phiw])
                    else:
                        phi = max(phic,phie)
            else:
                if x[0] < self.xsw[0]:
                    phi = phic
                else:
                    if x[0] < self.xse[0]:
                        phi = min(phisw,phise)
                    else:
                        phi = phic
            return -phi

analyticalSolution = {0:Zalesak2D(L=L,
                                  center=[0.0,
                                          0.5],
                                  radius=0.25,
                                  zalesak=True)}

class UnitSquareRotation(NCLS.Coefficients):
    from proteus.ctransportCoefficients import unitSquareRotationEvaluate
    from proteus.ctransportCoefficients import unitSquareRotationLevelSetEvaluate
    def __init__(self,useHJ=False,epsFact=1.5,checkMass=False,
                 RD_model=None,
                 useMetrics=0.0,sc_uref=1.0,sc_beta=1.0):
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
        NCLS.Coefficients.__init__(self)
        self.checkMass=checkMass
        self.useMetrics = 0.0
	self.sc_uref=1.0
	self.sc_beta=1.0
    def attachModels(self,modelList):
        self.model = modelList[0]
	self.u_old_dof = numpy.copy(self.model.u[0].dof)
        self.q_v = numpy.zeros(self.model.q[('dH',0,0)].shape,'d')
        self.ebqe_v = numpy.zeros(self.model.ebqe[('dH',0,0)].shape,'d')
        self.unitSquareRotationLevelSetEvaluate(self.model.timeIntegration.tLast,
                                              self.model.q['x'],
                                              self.model.q[('u',0)],self.model.q[('grad(u)',0)],
                                              self.model.q[('m',0)],self.model.q[('dm',0,0)],
                                              self.model.q[('dH',0,0)],self.model.q[('dH',0,0)],
                                              self.model.q[('H',0)],self.q_v)
        self.model.q[('velocity',0)]=self.q_v
        self.model.ebqe[('velocity',0)]=self.ebqe_v
        if self.RD_modelIndex != None:
            self.rdModel = modelList[self.RD_modelIndex]
        else:
            self.rdModel = self.model
    def preStep(self,t,firstStep=False):
        self.q_v[...,0]  = -2.0*math.pi*self.model.q['x'][...,1]
        self.q_v[...,1]  =  2.0*math.pi*self.model.q['x'][...,0]
        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
       	self.u_old_dof = numpy.copy(self.model.u[0].dof)
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        pass

if applyRedistancing:
    RD_model=1
else:
    RD_model=None

coefficients = UnitSquareRotation(useHJ=True,epsFact=epsFactHeaviside,checkMass=checkMass,RD_model=RD_model,useMetrics=useMetrics)

coefficients.variableNames=['u']

def getDBC(x,flag):
    pass

def zeroInflow(x):
    return lambda x,t: 0.0

dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
