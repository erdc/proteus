from proteus import *
from proteus.default_p import *
from math import *
from rotation2D import *
from proteus.mprans import NCLS
#import Profiling

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name = soname + "_ls"

nd = 2


class OscillatingRotation2D:
    # cek changed to put sphere inside arbitrary box with dimensions in L
    def __init__(self):
        self.radius = 0.1 * min(width_and_hight)
        self.xc = lower_left_cornor[0] + 0.5 * width_and_hight[0]
        self.yc = lower_left_cornor[1] + 0.6 * width_and_hight[1]

    def uOfXT(self, x, t):
        return 1.0 - numpy.tanh(((x[0] - self.xc)**2 + (x[1] - self.yc)**2) / self.radius / self.radius - 1)


class OscillatingRotation2Dcylinder:
    # cek changed to put sphere inside arbitrary box with dimensions in L
    def __init__(self, L):
        self.radius = 0.15 * L[0]
        self.xc = 0.5 * L[0]
        self.yc = 0.75 * L[1]

    def uOfXT(self, x, t):
        return self.radius - math.sqrt((x[0] - self.xc)**2 + (x[1] - self.yc)**2)


analyticalSolution = {0: OscillatingRotation2D()}


class UnitSquareRotation(NCLS.Coefficients):

    from proteus.ctransportCoefficients import ncLevelSetCoefficientsEvaluate

    def __init__(self, useHJ=True, epsFact=1.5, checkMass=False,
                 RD_model=None,
                 useMetrics=0.0, sc_uref=1.0, sc_beta=1.0):
        self.waterline_interval = -1
        self.epsFact = epsFact
        self.useHJ = useHJ
        self.RD_modelIndex = RD_model
        self.sc_uref = sc_uref
        self.sc_beta = sc_beta
        self.useMetrics = useMetrics
        mass = {0: {0: 'linear'}}
        advection = {0: {0: 'linear'}}
        diffusion = {}
        potential = {}
        reaction = {}
        if self.useHJ:
            hamiltonian = {0: {0: 'linear'}}
        else:
            hamiltonian = {}

        NCLS.Coefficients.__init__(self)

        self.checkMass = checkMass
        self.useMetrics = 0.0
        self.sc_uref = 1.0
        self.sc_beta = 1.0
        self.STABILIZATION_TYPE = ct.stablization

    def attachModels(self, modelList):
        self.model = modelList[0]
        self.u_old_dof = numpy.copy(self.model.u[0].dof)

        self.q_v = numpy.zeros(self.model.q[('dH', 0, 0)].shape, 'd')
        self.ebqe_v = numpy.zeros(self.model.ebqe[('dH', 0, 0)].shape, 'd')

        self.q_v[..., 0] = -2.0 * math.pi * \
            (self.model.q['x'][..., 1] - rotation_center[0])
        self.q_v[..., 1] = 2.0 * math.pi * \
            (self.model.q['x'][..., 0] - rotation_center[1])

        self.model.q[('velocity', 0)] = self.q_v
        self.model.ebqe[('velocity', 0)] = self.ebqe_v

        if self.RD_modelIndex != None:
            # print self.RD_modelIndex,len(modelList)
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
    def preStep(self, t, firstStep=False):

        self.q_v[..., 0] = -2.0 * math.pi * \
            (self.model.q['x'][..., 1] - rotation_center[0])
        self.q_v[..., 1] = 2.0 * math.pi * \
            (self.model.q['x'][..., 0] - rotation_center[1])

        copyInstructions = {}
        return copyInstructions
#

    def postStep(self, t, firstStep=False):
        self.u_old_dof = numpy.copy(self.model.u[0].dof)

        self.model.q['dV_last'][:] = self.model.q['dV']

        copyInstructions = {}
        return copyInstructions
#

    def evaluate(self, t, c):
        pass

    def calculateResidual(self, *args):
        import method_gp as M
        M.getResidual(*args)


if applyRedistancing:
    RD_model = 1
else:
    RD_model = None
coefficients = UnitSquareRotation(epsFact=epsFactHeaviside,
                                  checkMass=checkMass, RD_model=RD_model, useMetrics=useMetrics)

coefficients.variableNames = ['u']

# now define the Dirichlet boundary conditions


def getDBC(x, flag):
    None
    # if (x[1] == 0.0):
    #    return lambda x,t: 0.0
    # if (x[0] == 0.0 or
    #    x[0] == 1.0 or
    #    x[1] == 0.0 or
    #    x[1] == 1.0):
    #    return lambda x,t: 0.0


def zeroInflow(x):
    return lambda x, t: 0.0
    # if (x[0] == 0.0 and x[1] <= 0.5):
    #     return lambda x,t: 0.0
    # if (x[0] == 1.0 and x[1] >= 0.5):
    #     return lambda x,t: 0.0
    # if (x[1] == 0.0 and x[0] >= 0.5):
    #     return lambda x,t: 0.0
    # if (x[1] == 1.0 and x[0] <= 0.5):
    #     return lambda x,t: 0.0


dirichletConditions = {0: getDBC}

initialConditions = {0: analyticalSolution[0]}

fluxBoundaryConditions = {0: 'outFlow'}


def zeroadv(x):
    return lambda x, t: 0.0


advectiveFluxBoundaryConditions = {}
#advectiveFluxBoundaryConditions =  {0:zeroadv}


diffusiveFluxBoundaryConditions = {0: {}}

# @}
