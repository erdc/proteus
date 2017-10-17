from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from burgers2D import *
from proteus.mprans import VOF

name = soname + "_vof"

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

LevelModelType = VOF.LevelModel

# \ingroup test
#\file vof_vortex_2d_p.py
#
# \todo finish vof_vortex_2d_p.py

# if applyRedistancing:
#     coefficients = VOF.Coefficients(LS_model=0, V_model=0, RD_model=1, ME_model=2, checkMass=checkMass,
#                                     epsFact=epsFact_vof, useMetrics=useMetrics)
# elif not onlyVOF:
#     coefficients = VOF.Coefficients(LS_model=0, V_model=0, RD_model=None, ME_model=1, checkMass=checkMass,
#                                     epsFact=epsFact_vof, useMetrics=useMetrics)
# else:


# coefficients = VOF.Coefficients(RD_model=None, ME_model=0, checkMass=checkMass,
#                                 epsFact=epsFact_vof, useMetrics=useMetrics)
class UnitSquareRotationCoeff(VOF.Coefficients):

    def __init__(self, epsFact=1.5, checkMass=False,
                 RD_model=None,
                 useMetrics=0.0):

        self.waterline_interval = -1
        self.epsFact = epsFact
        self.ME_model = 0
        self.RD_modelIndex = RD_model
        self.useMetrics = useMetrics

        mass = {0: {0: 'linear'}}
        advection = {0: {0: 'linear'}}
        diffusion = {}
        potential = {}
        reaction = {}

        hamiltonian = {}

        VOF.Coefficients.__init__(self)

        self.checkMass = checkMass
        self.useMetrics = 0.0
        self.sc_uref = 1.0
        self.sc_beta = 1.0

        self.STABILIZATION_TYPE = ct.stablization
        self.FCT = (ct.stablization == 2)

    def attachModels(self, modelList):
        self.model = modelList[0]
        self.u_old_dof = numpy.copy(self.model.u[0].dof)

        self.q_v = numpy.zeros(self.model.q[('grad(u)', 0)].shape, 'd')
        self.ebqe_v = numpy.zeros(self.model.ebqe[('grad(u)', 0)].shape, 'd')

        self.ebqe_phi = numpy.zeros([1], 'd')
        self.epsFact = 0

        self.q_v[:, :, 0] = 0.5 * self.model.q[('u', 0)]
        self.q_v[:, :, 1] = 0.5 * self.model.q[('u', 0)]
        self.ebqe_v[:, :, 0] = 0.5 * self.model.ebqe[('u', 0)]
        self.ebqe_v[:, :, 1] = 0.5 * self.model.ebqe[('u', 0)]

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
        self.q_v[:, :, 0] = 0.5 * self.model.q[('u', 0)]
        self.q_v[:, :, 1] = 0.5 * self.model.q[('u', 0)]
        self.ebqe_v[:, :, 0] = 0.5 * self.model.ebqe[('u', 0)]
        self.ebqe_v[:, :, 1] = 0.5 * self.model.ebqe[('u', 0)]

        copyInstructions = {}
        return copyInstructions

    def postStep(self, t, firstStep=False):
        self.u_old_dof = numpy.copy(self.model.u[0].dof)
        copyInstructions = {}
        return copyInstructions

    def evaluate(self, t, c):
        pass


if applyRedistancing:
    RD_model = 1
else:
    RD_model = None
coefficients = UnitSquareRotationCoeff(epsFact=epsFactHeaviside,
                                       checkMass=checkMass, RD_model=RD_model, useMetrics=useMetrics)

coefficients.variableNames = ['u']


class Burgers2D:

    def true_value(self, x, y, t):
        u = 0
        if y > x:
            x, y = y, x
        a = x - y
        a0 = 1 - 0.5 * t
        if a <= a0:
            if 0 <= y < t:
                u = y / (t + 1e-10)
            elif t <= y < 0.5 * t + 1 - a:
                u = 1
            else:
                u = 0
        elif a <= 1:
            if 0 <= y < math.sqrt(2 * t * (1 - a)):
                u = y / (t + 1e-10)
            else:
                u = 0
        else:
            u = 0

        return u

    def uOfXT(self, XY, t):

        try:
            L = len(XY[0])
            u = numpy.zeros((L,), 'd')
            for i in xrange(L):
                u[i] = self.true_value(XY[0][i], XY[0][i])

        # when XY is a tuple of 3 numbers and XY[0] is a number
        except TypeError:
            u = self.true_value(XY[0], XY[1], t)

        return u


analyticalSolution = {0: Burgers2D()}


def getDBC(x, flag):
    pass


def zeroInflow(x):
    return lambda x, t: 0.0


dirichletConditions = {0: getDBC}


class initial_condition:
    def uOfXT(self, x, t=0):
        return (x[0] >= 0) * (x[0] <= 1) * (x[1] >= 0) * (x[1] <= 1)
# return ((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5) <=
# 0.25)


initialConditions = {0: initial_condition()}

fluxBoundaryConditions = {0: 'outFlow'}

# cek made no flux since v.n = 0 for this v


def getAFBC(x, flag):
    return lambda x, t: 0.0


advectiveFluxBoundaryConditions = {0: getAFBC}

diffusiveFluxBoundaryConditions = {0: {}}
