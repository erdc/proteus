from __future__ import absolute_import
from builtins import object
from proteus import *
from proteus.default_p import *
from math import *
try:
    from .rotation2D import *
except:
    from rotation2D import *
from proteus.mprans import NCLS
#import Profiling

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

## \page Tests Test Problems
# \ref ls_rotation_2d_p.py "Linear advection of a circular level set function in an oscillating rotation velocity field"
#

##\ingroup test
# \file la_rotation_2d_p.py
# @{
#  \brief Conservative linear advection of a circle signed distance function
#  in a oscillating rotation velocity field.
#
# \f{eqnarray*}
# \phi_t + \nabla \cdot (\vec u \phi) &=& 0 \\
# \Omega &=& [0,1] \times [0,1] \\
#  u^{x} &=& \cos(\pi t/8)\sin(2\pi y)\sin^2(\pi x) \\
#  u^{y} &=& -\cos(\pi t/8)\sin(2\pi x)\sin^{2}(\pi y) \\
#  \phi^{0} &=& \left(x-\frac{1}{2}\right)^2 + \left(y-\frac{3}{4}\right)^2 - 0.15^2
# \f}
# The solution should return to the initial condition at \f$T=8\f$.
# Outflow boundaries are applied on \f$\partial \Omega\f$.
#
#
# \image html  save_la_rotation_2d_dgp2_exact.jpg "exact solution, T=8.0"
# \image latex save_la_rotation_2d_dgp2_exact.eps "exact solution, T=8.0"
# \image html  save_la_rotation_2d_dgp2_phi.jpg "RKDG P^2 solution, Cr=0.1, L^2 error= 7.84e-3"
# \image latex save_la_rotation_2d_dgp2_phi.eps "RKDG $P^2$ solution, Cr=0.1, $L^2$ error= 7.84e-3"
#

class OscillatingRotation2D(object):
    #cek changed to put sphere inside arbitrary box with dimensions in L
    def __init__(self,L):
        self.radius = 0.25
        self.xc=0.0
        self.yc=0.5
    def uOfXT(self,x,t):
        return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
class OscillatingRotation2Dcylinder(object):
    #cek changed to put sphere inside arbitrary box with dimensions in L
    def __init__(self,L):
        self.radius = 0.15*L[0]
        self.xc=0.5*L[0]
        self.yc=0.75*L[1]
    def uOfXT(self,x,t):
        return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)

analyticalSolution = {0:OscillatingRotation2D(L)}

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
            #print self.RD_modelIndex,len(modelList)
            self.rdModel = modelList[self.RD_modelIndex]
        else:
            self.rdModel = self.model
        # if self.checkMass:
        #     self.m_pre = Norms.scalarSmoothedHeavisideDomainIntegral(self.epsFact,
        #                                                              self.model.mesh.elementDiametersArray,
        #                                                              self.model.q['dV'],
        #                                                              self.model.q[('m',0)],
        #                                                              self.model.mesh.nElements_owned)

        #     logEvent("Attach Models UnitSquareRotation: Phase  0 mass before NCLS step = %12.5e" % (self.m_pre,),level=2)
        #     self.totalFluxGlobal=0.0
        #     self.lsGlobalMassArray = [self.m_pre]
        #     self.lsGlobalMassErrorArray = [0.0]
        #     self.fluxArray = [0.0]
        #     self.timeArray = [self.model.timeIntegration.t]
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
