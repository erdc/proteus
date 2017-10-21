from proteus import *
from proteus.default_p import *
import math
from rotation2D import *
from proteus.mprans import NCLS
#import Profiling

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name = soname + "_ls"

nd = 2

#


class ProblemRotationALE2D:
    def __init__(self):
        self.radius = 0.1 * min(width_and_hight)
        self.xc = lower_left_cornor[0] + 0.5 * width_and_hight[0]
        self.yc = lower_left_cornor[1] + 0.7 * width_and_hight[1]
        self.z0 = rotation_center[0] + 1j * rotation_center[1]

    def uOfXT(self, x, t=0):
        # initial and final solution t=0 or 1
        return 1.0 - numpy.tanh(((x[0] - self.xc)**2 + (x[1] - self.yc)**2) / self.radius / self.radius - 1)

    @staticmethod
    def advection_function(x, y, t=0):
        return numpy.array([numpy.sin(numpy.pi * x) * numpy.cos(numpy.pi * y) * numpy.cos(2 * numpy.pi * t),
                            -numpy.cos(numpy.pi * x) * numpy.sin(numpy.pi * y) * numpy.cos(2 * numpy.pi * t)])

    @staticmethod
    def mesh_velocity(x, y, t=0):
        return np.stack((numpy.sin(numpy.pi * x) * numpy.cos(numpy.pi * y) * numpy.cos(2 * numpy.pi * t),
                         -numpy.cos(numpy.pi * x) * numpy.sin(numpy.pi *
                                                              y) * numpy.cos(2 * numpy.pi * t),
                         numpy.zeros_like(x, 'd')),
                        axis=1)


analyticalSolution = {0: ProblemRotationALE2D()}


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

    def save_dirichlet_dofs(self):
        mesh = self.model.mesh
        fes = self.model.u[0].femSpace
        self.dirichlet_bc_dofs = {'dof': [],
                                  'xyz': [], 'label': [], 'value': []}

        for eN in range(mesh.nElements_global):
            for k in range(fes.referenceFiniteElement.interpolationConditions.nQuadraturePoints):
                i = fes.referenceFiniteElement.interpolationConditions.quadrature2DOF_element(
                    k)
                dofN = fes.dofMap.l2g[eN, i]
                x = fes.interpolationPoints[eN, k]
                for ebN_element in range(mesh.nElementBoundaries_element):
                    if fes.referenceFiniteElement.interpolationConditions.definedOnLocalElementBoundary(k, ebN_element) == True:
                        ebN = mesh.elementBoundariesArray[eN, ebN_element]
                        materialFlag = mesh.elementBoundaryMaterialTypes[ebN]
                        if materialFlag in [1, 2, 3, 4]:
                            self.dirichlet_bc_dofs['dof'].append(dofN)
                            self.dirichlet_bc_dofs['xyz'].append(x)
                            self.dirichlet_bc_dofs['label'].append(
                                materialFlag)
                            self.dirichlet_bc_dofs['value'].append(0)

        self.dirichlet_bc_dofs['dof'] = numpy.asarray(
            self.dirichlet_bc_dofs['dof'], 'i')
        self.dirichlet_bc_dofs['xyz'] = numpy.asarray(
            self.dirichlet_bc_dofs['xyz'], 'd')
        self.dirichlet_bc_dofs['label'] = numpy.asarray(
            self.dirichlet_bc_dofs['label'], 'd')
        self.dirichlet_bc_dofs['value'] = numpy.asarray(
            self.dirichlet_bc_dofs['value'], 'd')

    def attachModels(self, modelList):
        self.model = modelList[0]
        self.mesh = self.model.mesh
        # It must be related to self.model.u_dof_old so that postStep
        # can upate it corectly.
        self.u_old_dof = numpy.copy(self.model.u[0].dof)

        self.q_v = numpy.zeros(self.model.q[('dH', 0, 0)].shape, 'd')
        self.ebqe_v = numpy.zeros(self.model.ebqe[('dH', 0, 0)].shape, 'd')

#         self.q_v[..., 0] = -2.0 * math.pi * \
#             (self.model.q['x'][..., 1] - rotation_center[0])
#         self.q_v[..., 1] = 2.0 * math.pi * \
#             (self.model.q['x'][..., 0] - rotation_center[1])
        self.q_v[..., 0] = 1.0
        self.q_v[..., 1] = 0.0

        self.model.q[('velocity', 0)] = self.q_v
        self.model.ebqe[('velocity', 0)] = self.ebqe_v

        if self.RD_modelIndex != None:
            # print self.RD_modelIndex,len(modelList)
            self.rdModel = modelList[self.RD_modelIndex]
        else:
            self.rdModel = self.model

        self.save_dirichlet_dofs()
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

        #         self.q_v[..., 0] = -2.0 * math.pi * \
        #             (self.model.q['x'][..., 1] - rotation_center[0])
        #         self.q_v[..., 1] = 2.0 * math.pi * \
        #             (self.model.q['x'][..., 0] - rotation_center[1])
        self.q_v[..., 0] = 1.0
        self.q_v[..., 1] = 0.0

        self.mesh.nodeVelocityArray[:] = analyticalSolution[0].mesh_velocity(
            self.mesh.nodeArray[:, 0], self.mesh.nodeArray[:, 1], t)

        copyInstructions = {}
        return copyInstructions
#

    def postStep(self, t, firstStep=False):
        #         import pdb
        #         pdb.set_trace()
        # This is called from proteus/NumericalSolution.py(1524)postStep()
        # Serious error it should be [:]
        self.u_old_dof[:] = numpy.copy(self.model.u[0].dof)

        self.model.q['dV_last'][:] = self.model.q['dV']

        # It should be model.stepController.dt_model
        # But now there is only 1 substep
        self.mesh.nodeArray[:] += self.model.timeIntegration.dt * \
            self.mesh.nodeVelocityArray
        print ">>>>>>>>>>>>>>>>>.moving mesh:", self.model.timeIntegration.dt
        copyInstructions = {}
        return copyInstructions
#

    def evaluate(self, t, c):
        pass

    def calculateResidual(self, *args):
        #         import pdb
        #         pdb.set_trace()
        #import method_LxF as M
        args[49][:] = args[33]
        args[43][:] = 40
        args[44][:] = 50
        return None

        import method_GP_ALE as M
        M.getResidual(
            *args, ad_function=analyticalSolution[0].advection_function)

        args[49][self.dirichlet_bc_dofs['dof']
                 ] = self.dirichlet_bc_dofs['value']


if applyRedistancing:
    RD_model = 1
else:
    RD_model = None
coefficients = UnitSquareRotation(epsFact=epsFactHeaviside,
                                  checkMass=checkMass, RD_model=RD_model, useMetrics=useMetrics)

coefficients.variableNames = ['u']

# now define the Dirichlet boundary conditions


def getDBC(x, flag):
    if (x[0] == 0.0 or
        x[0] == 1.0 or
        x[1] == 0.0 or
            x[1] == 1.0):
        return lambda x, t: 0.0


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
