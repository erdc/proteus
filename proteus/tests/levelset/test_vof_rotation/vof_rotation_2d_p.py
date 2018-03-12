from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from rotation2D import *
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

    def __init__(self, epsFact=1.5, ME_model=0, checkMass=False,
                 RD_model=None,
                 useMetrics=0.0):

        self.waterline_interval = -1
        self.epsFact = epsFact
        self.ME_model = ME_model
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


def Heaviside(phi):
    if phi > 0:
        return 1.0
    elif phi < 0:
        return 0.0
    else:
        return 0.5


class Vortex_phi:
    def __init__(self):
        self.radius = 0.1 * min(width_and_hight)
        self.center = [lower_left_cornor[0] + 0.5 * width_and_hight[0],
                       lower_left_cornor[1] + 0.6 * width_and_hight[1]]

    def uOfX(self, X):
        dx = X[0] - self.center[0]
        dy = X[1] - self.center[1]
        dBubble = self.radius - sqrt(dx**2 + dy**2)
        # Heaviside(dBubble)
        return smoothedHeaviside(epsFactHeaviside * he, dBubble)
    # end

    def uOfXT(self, X, t):
        return self.uOfX(X)
    # end
# end Vortex_phi


class Vortex_phi_cylinder:
    def __init__(self, center=[0.5, 0.75, 0.5], radius=0.15):
        self.radius = radius
        self.center = center

    def uOfX(self, X):
        dx = X[0] - self.center[0]
        dy = X[1] - self.center[1]
        dBubble = self.radius - sqrt(dx**2 + dy**2)
        # Heaviside(dBubble)
        return smoothedHeaviside(epsFactHeaviside * he, dBubble)
    # end

    def uOfXT(self, X, t):
        return self.uOfX(X)
    # end
# end Vortex_phi


def getDBC(x, flag):
    pass


dirichletConditions = {0: getDBC}

analyticalSolution = {0: Vortex_phi()}
initialConditions = {0: Vortex_phi()}

fluxBoundaryConditions = {0: 'outFlow'}

# cek made no flux since v.n = 0 for this v


def getAFBC(x, flag):
    return lambda x, t: 0.0


advectiveFluxBoundaryConditions = {0: getAFBC}

diffusiveFluxBoundaryConditions = {0: {}}
