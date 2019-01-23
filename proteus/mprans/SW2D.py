from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import proteus
from proteus.mprans.cSW2D import *
from proteus.Comm import globalMax

class SubgridError(proteus.SubgridError.SGE_base):
    def __init__(self, coefficients, nd, lag=False, nStepsToDelay=0, hFactor=1.0):
        proteus.SubgridError.SGE_base.__init__(self, coefficients, nd, lag)
        self.hFactor = hFactor
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0
        if self.lag:
            logEvent("SW2D.SubgridError: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay = 1
            self.lag = False

    def initializeElementQuadrature(self, mesh, t, cq):
        import copy
        self.cq = cq
        self.v_last = self.cq[('velocity', 0)]

    def updateSubgridErrorHistory(self, initializationPhase=False):
        self.nSteps += 1
        if self.lag:
            self.v_last[:] = self.cq[('velocity', 0)]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            logEvent("SW2D.SubgridError: switched to lagged subgrid error")
            self.lag = True
            self.v_last = self.cq[('velocity', 0)].copy()

    def calculateSubgridError(self, q):
        pass


class NumericalFlux(proteus.NumericalFlux.ShallowWater_2D):
    hasInterior = False

    def __init__(self, vt, getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 h_eps=1.0e-8,
                 tol_u=1.0e-8):
        proteus.NumericalFlux.ShallowWater_2D.__init__(self, vt, getPointwiseBoundaryConditions,
                                                       getAdvectiveFluxBoundaryConditions,
                                                       getDiffusiveFluxBoundaryConditions,
                                                       getPeriodicBoundaryConditions,
                                                       h_eps,
                                                       tol_u)
        self.penalty_constant = 2.0
        self.includeBoundaryAdjoint = True
        self.boundaryAdjoint_sigma = 1.0
        self.hasInterior = False


class ShockCapturing(proteus.ShockCapturing.ShockCapturing_base):
    def __init__(self, coefficients, nd, shockCapturingFactor=0.25, lag=False, nStepsToDelay=3):
        proteus.ShockCapturing.ShockCapturing_base.__init__(self, coefficients, nd, shockCapturingFactor, lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps = 0
        if self.lag:
            logEvent("SW2D.ShockCapturing: lagging requested but must lag the first step; switching lagging off and delaying")
            self.nStepsToDelay = 1
            self.lag = False

    def initializeElementQuadrature(self, mesh, t, cq):
        self.mesh = mesh
        self.numDiff = {}
        self.numDiff_last = {}
        for ci in range(3):
            self.numDiff[ci] = cq[('numDiff', ci, ci)]
            self.numDiff_last[ci] = cq[('numDiff', ci, ci)]

    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(3):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay is not None and self.nSteps > self.nStepsToDelay:
            logEvent("SW2D.ShockCapturing: switched to lagged shock capturing")
            self.lag = True
            for ci in range(3):
                self.numDiff_last[ci] = self.numDiff[ci].copy()
        logEvent("SW2D: max numDiff_0 %e numDiff_1 %e numDiff_2 %e" % (globalMax(self.numDiff_last[0].max()),
                                                                       globalMax(self.numDiff_last[1].max()),
                                                                       globalMax(self.numDiff_last[2].max())))


class Coefficients(proteus.TransportCoefficients.TC_base):
    """
    The coefficients for the shallow water equations
    """

    def __init__(self,
                 nu=1.004e-6,
                 g=9.8,
                 nd=2,
                 sd=True,
                 movingDomain=False,
                 useRBLES=0.0,
                 useMetrics=0.0):
        self.useRBLES = useRBLES
        self.useMetrics = useMetrics
        self.sd = sd
        self.nu = nu
        self.g = g
        self.nd = nd
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction = {}
        hamiltonian = {}
        if nd == 2:
            variableNames = ['h', 'u', 'v']
            mass = {0: {0: 'linear'},
                    1: {0: 'linear', 1: 'linear'},
                    2: {0: 'linear', 2: 'linear'}}
            advection = {0: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'},
                         1: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'},
                         2: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear'}}
            diffusion = {1: {1: {1: 'constant'}, 2: {2: 'constant'}},
                         2: {2: {2: 'constant'}, 1: {1: 'constant'}}}
            sdInfo = {(1, 1): (numpy.array([0, 1, 2], dtype='i'),
                               numpy.array([0, 1], dtype='i')),
                      (1, 2): (numpy.array([0, 0, 1], dtype='i'),
                               numpy.array([0], dtype='i')),
                      (2, 2): (numpy.array([0, 1, 2], dtype='i'),
                               numpy.array([0, 1], dtype='i')),
                      (2, 1): (numpy.array([0, 1, 1], dtype='i'),
                               numpy.array([1], dtype='i'))}
            potential = {1: {1: 'u'},
                         2: {2: 'u'}}
            reaction = {1: {0: 'linear'},
                        2: {0: 'linear'}}
            TC_base.__init__(self,
                             3,
                             mass,
                             advection,
                             diffusion,
                             potential,
                             reaction,
                             hamiltonian,
                             variableNames,
                             sparseDiffusionTensors=sdInfo,
                             useSparseDiffusion=sd,
                             movingDomain=movingDomain)
            self.vectorComponents = [1, 2]

    def attachModels(self, modelList):
        pass

    def initializeMesh(self, mesh):
        pass

    def initializeElementQuadrature(self, t, cq):
        pass

    def initializeElementBoundaryQuadrature(self, t, cebq, cebq_global):
        pass

    def initializeGlobalExteriorElementBoundaryQuadrature(self, t, cebqe):
        pass

    def updateToMovingDomain(self, t, c):
        pass

    def evaluate(self, t, c):
        pass


class LevelModel(proteus.Transport.OneLevelTransport):
    nCalls = 0

    def __init__(self,
                 uDict,
                 phiDict,
                 testSpaceDict,
                 matType,
                 dofBoundaryConditionsDict,
                 dofBoundaryConditionsSetterDict,
                 coefficients,
                 elementQuadrature,
                 elementBoundaryQuadrature,
                 fluxBoundaryConditionsDict=None,
                 advectiveFluxBoundaryConditionsSetterDict=None,
                 diffusiveFluxBoundaryConditionsSetterDictDict=None,
                 stressTraceBoundaryConditionsSetterDictDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='SW2D',
                 reuse_trial_and_test_quadrature=True,
                 sd=True,
                 movingDomain=False):
        self.postProcessing = False  # this is a hack to test the effect of post-processing
        #
        # set the objects describing the method and boundary conditions
        #
        self.movingDomain = movingDomain
        self.tLast_mesh = None
        #
        # cek todo clean up these flags in the optimized version
        self.bcsTimeDependent = options.bcsTimeDependent
        self.bcsSet = False
        self.name = name
        self.sd = sd
        self.lowmem = True
        self.timeTerm = True  # allow turning off  the  time derivative
        self.testIsTrial = True
        self.phiTrialIsTrial = True
        self.u = uDict
        self.Hess = False
        if isinstance(self.u[0].femSpace, C0_AffineQuadraticOnSimplexWithNodalBasis):
            self.Hess = True
        self.ua = {}  # analytical solutions
        self.phi = phiDict
        self.dphi = {}
        self.matType = matType
        # mwf try to reuse test and trial information across components if spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature  # True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1, coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        # Simplicial Mesh
        self.mesh = self.u[0].femSpace.mesh  # assume the same mesh for  all components for now
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        self.dirichletNodeSetList = None  # explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.coefficients = coefficients
        # cek hack? give coefficients a bathymetriy array
        import copy
        self.coefficients.b = self.u[0].copy()
        self.coefficients.b.name = 'b'
        self.coefficients.b.dof.fill(0.0)
        #
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        self.conservativeFlux = conservativeFluxDict  # no velocity post-processing for now
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        # determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        # cek come back
        if self.stabilization is not None:
            for ci in range(self.nc):
                if ci in coefficients.mass:
                    for flag in list(coefficients.mass[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.advection:
                    for flag in list(coefficients.advection[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.diffusion:
                    for diffusionDict in list(coefficients.diffusion[ci].values()):
                        for flag in list(diffusionDict.values()):
                            if flag != 'constant':
                                self.stabilizationIsNonlinear = True
                if ci in coefficients.potential:
                    for flag in list(coefficients.potential[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.reaction:
                    for flag in list(coefficients.reaction[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if ci in coefficients.hamiltonian:
                    for flag in list(coefficients.hamiltonian[ci].values()):
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
        # determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci in range(self.nc):
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux is not None) or
                                                 (numericalFluxType is not None) or
                                                 (self.fluxBoundaryConditions[ci] == 'outFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'mixedFlow') or
                                                 (self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        # calculate some dimensions
        #
        self.nSpace_global = self.u[0].femSpace.nSpace_global  # assume same space dim for all variables
        self.nDOF_trial_element = [u_j.femSpace.max_nDOF_element for u_j in list(self.u.values())]
        self.nDOF_phi_trial_element = [phi_k.femSpace.max_nDOF_element for phi_k in list(self.phi.values())]
        self.n_phi_ip_element = [phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for phi_k in list(self.phi.values())]
        self.nDOF_test_element = [femSpace.max_nDOF_element for femSpace in list(self.testSpace.values())]
        self.nFreeDOF_global = [dc.nFreeDOF_global for dc in list(self.dirichletConditions.values())]
        self.nVDOF_element = sum(self.nDOF_trial_element)
        self.nFreeVDOF_global = sum(self.nFreeDOF_global)
        #
        NonlinearEquation.__init__(self, self.nFreeVDOF_global)
        #
        # build the quadrature point dictionaries from the input (this
        # is just for convenience so that the input doesn't have to be
        # complete)
        #
        elementQuadratureDict = {}
        elemQuadIsDict = isinstance(elementQuadrature, dict)
        if elemQuadIsDict:  # set terms manually
            for I in self.coefficients.elementIntegralKeys:
                if I in elementQuadrature:
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization is not None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if I in elementQuadrature:
                        elementQuadratureDict[('stab',) + I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',) + I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',) + I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if ('numDiff', ci, ci) in elementQuadrature:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[('numDiff', ci, ci)]
                    else:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature
        if massLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('m', ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        if reactionLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('r', ci)] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        elementBoundaryQuadratureDict = {}
        if isinstance(elementBoundaryQuadrature, dict):  # set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if I in elementBoundaryQuadrature:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature['default']
        else:
            for I in self.coefficients.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureDict[I] = elementBoundaryQuadrature
        #
        # find the union of all element quadrature points and
        # build a quadrature rule for each integral that has a
        # weight at each point in the union
        # mwf include tag telling me which indices are which quadrature rule?
        (self.elementQuadraturePoints, self.elementQuadratureWeights,
         self.elementQuadratureRuleIndeces) = Quadrature.buildUnion(elementQuadratureDict)
        self.nQuadraturePoints_element = self.elementQuadraturePoints.shape[0]
        self.nQuadraturePoints_global = self.nQuadraturePoints_element * self.mesh.nElements_global
        #
        # Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints,
         self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[0]
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global *
                                                        self.mesh.nElementBoundaries_element *
                                                        self.nElementBoundaryQuadraturePoints_elementBoundary)

#        if isinstance(self.u[0].femSpace,C0_AffineLinearOnSimplexWithNodalBasis):
#            print self.nQuadraturePoints_element
#            if self.nSpace_global == 3:
#                assert(self.nQuadraturePoints_element == 5)
#            elif self.nSpace_global == 2:
#                assert(self.nQuadraturePoints_element == 6)
#            elif self.nSpace_global == 1:
#                assert(self.nQuadraturePoints_element == 3)
#
#            print self.nElementBoundaryQuadraturePoints_elementBoundary
#            if self.nSpace_global == 3:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 2:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
#            elif self.nSpace_global == 1:
#                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)
        #
        # simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q = {}
        self.ebq = {}
        self.ebq_global = {}
        self.ebqe = {}
        self.phi_ip = {}
        # mesh
        self.h_dof_sge = self.u[0].dof.copy()
        self.u_dof_sge = self.u[1].dof.copy()
        self.v_dof_sge = self.u[2].dof.copy()
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary, 3), 'd')
        self.ebq_global[('totalFlux', 0)] = numpy.zeros((self.mesh.nElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebq_global[('velocityAverage', 0)] = numpy.zeros((self.mesh.nElementBoundaries_global,
                                                               self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.q[('u', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m', 0)] = self.q[('u', 0)]
        self.q[('m', 1)] = self.q[('u', 1)]
        self.q[('m', 2)] = self.q[('u', 2)]
        self.q[('m_last', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('f', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('velocity', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('cfl', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 0, 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 1, 1)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('numDiff', 2, 2)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.ebqe[('u', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc_flag', 0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 1, 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 2, 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 1, 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('penalty')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 2, 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('velocity', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('velocity', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('velocity', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set([('u', ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature = set()
        self.tensors_elementBoundaryQuadrature = set()
        #
        # show quadrature
        #
        logEvent("Dumping quadrature shapes for model %s" % self.name, level=9)
        logEvent("Element quadrature array (q)", level=9)
        for (k, v) in list(self.q.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Element boundary quadrature (ebq)", level=9)
        for (k, v) in list(self.ebq.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Global element boundary quadrature (ebq_global)", level=9)
        for (k, v) in list(self.ebq_global.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Exterior element boundary quadrature (ebqe)", level=9)
        for (k, v) in list(self.ebqe.items()):
            logEvent(str((k, v.shape)), level=9)
        logEvent("Interpolation points for nonlinear diffusion potential (phi_ip)", level=9)
        for (k, v) in list(self.phi_ip.items()):
            logEvent(str((k, v.shape)), level=9)
        #
        # allocate residual and Jacobian storage
        #
        #
        # allocate residual and Jacobian storage
        #
        self.elementResidual = [numpy.zeros(
            (self.mesh.nElements_global,
             self.nDOF_test_element[ci]),
            'd')]
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
            self.inflowBoundaryBC[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,), 'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nDOF_trial_element[cj]), 'd')
            self.inflowFlux[cj] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        # identify the internal nodes this is ought to be in mesh
        # \todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global = self.mesh.elementBoundaryElementsArray[ebN, 0]
            ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[ebN, 0]
            for i in range(self.mesh.nNodes_element):
                if i != ebN_element:
                    I = self.mesh.elementNodesArray[eN_global, i]
                    self.internalNodes -= set([I])
        self.nNodes_internal = len(self.internalNodes)
        self.internalNodesArray = numpy.zeros((self.nNodes_internal,), 'i')
        for nI, n in enumerate(self.internalNodes):
            self.internalNodesArray[nI] = n
        #
        del self.internalNodes
        self.internalNodes = None
        logEvent("Updating local to global mappings", 2)
        self.updateLocal2Global()
        logEvent("Building time integration object", 2)
        logEvent(memory("inflowBC, internalNodes,updateLocal2Global", "OneLevelTransport"), level=4)
        # mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(self, integrateInterpolationPoints=True)
        else:
            self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        logEvent(memory("TimeIntegration", "OneLevelTransport"), level=4)
        logEvent("Calculating numerical quadrature formulas", 2)
        self.calculateQuadrature()
        self.setupFieldStrides()

        comm = Comm.get()
        self.comm = comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions, "You must use a numerical flux to apply weak boundary conditions for parallel runs"

        logEvent(memory("stride+offset", "OneLevelTransport"), level=4)
        if numericalFluxType is not None:
            if options is None or options.periodicDirichletConditions is None:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(self,
                                                       dofBoundaryConditionsSetterDict,
                                                       advectiveFluxBoundaryConditionsSetterDict,
                                                       diffusiveFluxBoundaryConditionsSetterDictDict,
                                                       options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        # set penalty terms
        # cek todo move into numerical flux initialization
        if 'penalty' in self.ebq_global:
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN, k] = old_div(self.numericalFlux.penalty_constant, \
                        (self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power))
        # penalty term
        # cek move  to Numerical flux initialization
        if 'penalty' in self.ebqe:
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE, k] = old_div(self.numericalFlux.penalty_constant, \
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        logEvent(memory("numericalFlux", "OneLevelTransport"), level=4)
        self.elementEffectiveDiametersArray = self.mesh.elementInnerDiametersArray
        # use post processing tools to get conservative fluxes, None by default
        if self.postProcessing:
            self.q[('v', 0)] = self.tmpvt.q[('v', 0)]
            self.ebq[('v', 0)] = self.tmpvt.ebq[('v', 0)]
            self.ebq[('w', 0)] = self.tmpvt.ebq[('w', 0)]
            self.ebq['sqrt(det(g))'] = self.tmpvt.ebq['sqrt(det(g))']
            self.ebq['n'] = self.tmpvt.ebq['n']
            self.ebq[('dS_u', 0)] = self.tmpvt.ebq[('dS_u', 0)]
            self.ebqe['dS'] = self.tmpvt.ebqe['dS']
            self.ebqe['n'] = self.tmpvt.ebqe['n']
            self.ebq_global['n'] = self.tmpvt.ebq_global['n']
            self.ebq_global['x'] = self.tmpvt.ebq_global['x']
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)
        logEvent(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        for ci, fbcObject in list(self.fluxBoundaryConditionsObjectsDict.items()):
            self.ebqe[('advectiveFlux_bc_flag', ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc', ci)].shape, 'i')
            for t, g in list(fbcObject.advectiveFluxBoundaryConditionsDict.items()):
                if ci in self.coefficients.advection:
                    self.ebqe[('advectiveFlux_bc', ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1
            for ck, diffusiveFluxBoundaryConditionsDict in list(fbcObject.diffusiveFluxBoundaryConditionsDictDict.items()):
                self.ebqe[('diffusiveFlux_bc_flag', ck, ci)] = numpy.zeros(self.ebqe[('diffusiveFlux_bc', ck, ci)].shape, 'i')
                for t, g in list(diffusiveFluxBoundaryConditionsDict.items()):
                    self.ebqe[('diffusiveFlux_bc', ck, ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag', ck, ci)][t[0], t[1]] = 1
        # self.numericalFlux.setDirichletValues(self.ebqe)
        if self.movingDomain:
            self.MOVING_DOMAIN = 1.0
        else:
            self.MOVING_DOMAIN = 0.0
        # cek hack
        self.movingDomain = False
        self.MOVING_DOMAIN = 0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(self.mesh.nodeArray.shape, 'd')
        # cek/ido todo replace python loops in modules with optimized code if possible/necessary
        self.forceStrongConditions = True  # False
        self.dirichletConditionsForceDOF = {}
        if self.forceStrongConditions:
            for cj in range(self.nc):
                self.dirichletConditionsForceDOF[cj] = DOFBoundaryConditions(
                    self.u[cj].femSpace, dofBoundaryConditionsSetterDict[cj], weakDirichletConditions=False)

        compKernelFlag = 0
        # if self.coefficients.useConstantH:
        #    self.elementDiameter = self.mesh.elementDiametersArray.copy()
        #    self.elementDiameter[:] = max(self.mesh.elementDiametersArray)
        # else:
        self.elementDiameter = self.mesh.elementDiametersArray
        print(self.nSpace_global, " nSpace_global")
        self.sw2d = cSW2D_base(self.nSpace_global,
                               self.nQuadraturePoints_element,
                               self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                               self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                               self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                               self.nElementBoundaryQuadraturePoints_elementBoundary,
                               compKernelFlag)
        try_supg = False
        if options is not None and 'try_supg_stabilization' in dir(options):
            try_supg = options.try_supg_stabilization
            logEvent("setting try_supg_stabilization from options= {0}".format(try_supg), level=1)
        if try_supg:
            self.calculateResidual = self.sw2d.calculateResidual_supg
            self.calculateJacobian = self.sw2d.calculateJacobian_supg
        else:
            self.calculateResidual = self.sw2d.calculateResidual
            self.calculateJacobian = self.sw2d.calculateJacobian

    def getResidual(self, u, r):
        """
        Calculate the element residuals and add in to the global residual
        """

        # Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        # cek todo put in logic to skip if BC's don't depend on t or u
        if self.bcsTimeDependent or not self.bcsSet:
            self.bcsSet = True
            # Dirichlet boundary conditions
            self.numericalFlux.setDirichletValues(self.ebqe)
            # Flux boundary conditions
            for ci, fbcObject in list(self.fluxBoundaryConditionsObjectsDict.items()):
                for t, g in list(fbcObject.advectiveFluxBoundaryConditionsDict.items()):
                    if ci in self.coefficients.advection:
                        self.ebqe[('advectiveFlux_bc', ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                        self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1
                for ck, diffusiveFluxBoundaryConditionsDict in list(fbcObject.diffusiveFluxBoundaryConditionsDictDict.items()):
                    for t, g in list(diffusiveFluxBoundaryConditionsDict.items()):
                        self.ebqe[('diffusiveFlux_bc', ck, ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                        self.ebqe[('diffusiveFlux_bc_flag', ck, ci)][t[0], t[1]] = 1
        r.fill(0.0)
        self.Ct_sge = 4.0
        self.Cd_sge = 144.0

        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.items()):
                    self.u[cj].dof[dofN] = g(self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN], self.timeIntegration.t)
        #import pdb
        # pdb.set_trace()

        self.calculateResidual(  # element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            # physics
            self.elementDiameter,  # mesh.elementDiametersArray,
            self.mesh.nElements_global,
            self.coefficients.useRBLES,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.coefficients.nu,
            self.coefficients.g,
            self.shockCapturing.shockCapturingFactor,
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.coefficients.b.dof,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.h_dof_sge,
            self.u_dof_sge,
            self.v_dof_sge,
            self.timeIntegration.m_tmp[0],
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.q[('f', 0)],
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.stabilization.v_last,
            self.q[('cfl', 0)],
            self.q[('numDiff', 0, 0)],
            self.q[('numDiff', 1, 1)],
            self.q[('numDiff', 2, 2)],
            self.shockCapturing.numDiff_last[0],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.coefficients.sdInfo[(1, 1)][0],
            self.coefficients.sdInfo[(1, 1)][1],
            self.coefficients.sdInfo[(1, 2)][0],
            self.coefficients.sdInfo[(1, 2)][1],
            self.coefficients.sdInfo[(2, 2)][0],
            self.coefficients.sdInfo[(2, 2)][1],
            self.coefficients.sdInfo[(2, 1)][0],
            self.coefficients.sdInfo[(2, 1)][1],
            self.offset[0],
            self.offset[1],
            self.offset[2],
            self.stride[0],
            self.stride[1],
            self.stride[2],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 1)],
            self.ebqe[('advectiveFlux_bc_flag', 2)],
            self.ebqe[('diffusiveFlux_bc_flag', 1, 1)],
            self.ebqe[('diffusiveFlux_bc_flag', 2, 2)],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.ebqe[('advectiveFlux_bc', 1)],
            self.ebqe[('advectiveFlux_bc', 2)],
            self.numericalFlux.ebqe[('u', 1)],
            self.ebqe[('diffusiveFlux_bc', 1, 1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u', 2)],
            self.ebqe[('diffusiveFlux_bc', 2, 2)],
            self.q[('velocity', 0)],
            self.ebqe[('velocity', 0)],
            self.ebq_global[('totalFlux', 0)],
            self.elementResidual[0])
        #import pdb
        # pdb.set_trace()
        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.items()):
                    r[self.offset[cj] + self.stride[cj] * dofN] = 0

        cflMax = globalMax(self.q[('cfl', 0)].max()) * self.timeIntegration.dt
        logEvent("Maximum CFL = " + str(cflMax), level=2)

        if self.stabilization:
            self.stabilization.accumulateSubgridMassHistory(self.q)
        logEvent("Global residual", level=9, data=r)
        # mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        self.calculateJacobian(  # element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.nodeVelocityArray,
            self.MOVING_DOMAIN,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            self.u[1].femSpace.psi,
            self.u[1].femSpace.grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[1].femSpace.psi_trace,
            self.u[1].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.elementDiameter,  # mesh.elementDiametersArray,
            self.mesh.nElements_global,
            self.coefficients.useRBLES,
            self.coefficients.useMetrics,
            self.timeIntegration.alpha_bdf,
            self.coefficients.nu,
            self.coefficients.g,
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.coefficients.b.dof,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.h_dof_sge,
            self.u_dof_sge,
            self.v_dof_sge,
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.stabilization.v_last,
            self.q[('cfl', 0)],
            self.shockCapturing.numDiff_last[0],
            self.shockCapturing.numDiff_last[1],
            self.shockCapturing.numDiff_last[2],
            self.coefficients.sdInfo[(1, 1)][0],
            self.coefficients.sdInfo[(1, 1)][1],
            self.coefficients.sdInfo[(1, 2)][0],
            self.coefficients.sdInfo[(1, 2)][1],
            self.coefficients.sdInfo[(2, 2)][0],
            self.coefficients.sdInfo[(2, 2)][1],
            self.coefficients.sdInfo[(2, 1)][0],
            self.coefficients.sdInfo[(2, 1)][1],
            self.csrRowIndeces[(0, 0)],
            self.csrColumnOffsets[(0, 0)],
            self.csrRowIndeces[(0, 1)],
            self.csrColumnOffsets[(0, 1)],
            self.csrRowIndeces[(0, 2)],
            self.csrColumnOffsets[(0, 2)],
            self.csrRowIndeces[(1, 0)],
            self.csrColumnOffsets[(1, 0)],
            self.csrRowIndeces[(1, 1)],
            self.csrColumnOffsets[(1, 1)],
            self.csrRowIndeces[(1, 2)],
            self.csrColumnOffsets[(1, 2)],
            self.csrRowIndeces[(2, 0)], self.csrColumnOffsets[(2, 0)],
            self.csrRowIndeces[(2, 1)], self.csrColumnOffsets[(2, 1)],
            self.csrRowIndeces[(2, 2)], self.csrColumnOffsets[(2, 2)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.isDOFBoundary[1],
            self.numericalFlux.isDOFBoundary[2],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 1)],
            self.ebqe[('advectiveFlux_bc_flag', 2)],
            self.ebqe[('diffusiveFlux_bc_flag', 1, 1)],
            self.ebqe[('diffusiveFlux_bc_flag', 2, 2)],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.ebqe[('advectiveFlux_bc', 1)],
            self.ebqe[('advectiveFlux_bc', 2)],
            self.numericalFlux.ebqe[('u', 1)],
            self.ebqe[('diffusiveFlux_bc', 1, 1)],
            self.ebqe[('penalty')],
            self.numericalFlux.ebqe[('u', 2)],
            self.ebqe[('diffusiveFlux_bc', 2, 2)],
            self.csrColumnOffsets_eb[(0, 0)],
            self.csrColumnOffsets_eb[(0, 1)],
            self.csrColumnOffsets_eb[(0, 2)],
            self.csrColumnOffsets_eb[(1, 0)],
            self.csrColumnOffsets_eb[(1, 1)],
            self.csrColumnOffsets_eb[(1, 2)],
            self.csrColumnOffsets_eb[(2, 0)],
            self.csrColumnOffsets_eb[(2, 1)],
            self.csrColumnOffsets_eb[(2, 2)])

        # Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            scaling = 1.0  # probably want to add some scaling to match non-dirichlet diagonals in linear system
            for cj in range(self.nc):
                for dofN in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.keys()):
                    global_dofN = self.offset[cj] + self.stride[cj] * dofN
                    for i in range(self.rowptr[global_dofN], self.rowptr[global_dofN + 1]):
                        if (self.colind[i] == global_dofN):
                            # print "RBLES forcing residual cj = %s dofN= %s global_dofN= %s was self.nzval[i]= %s now =%s " % (cj,dofN,global_dofN,self.nzval[i],scaling)
                            self.nzval[i] = scaling
                        else:
                            self.nzval[i] = 0.0
                            # print "RBLES zeroing residual cj = %s dofN= %s global_dofN= %s " % (cj,dofN,global_dofN)
        logEvent("Jacobian ", level=10, data=jacobian)
        # mwf decide if this is reasonable for solver statistics
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian

    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.tmpvt.calculateElementQuadrature()
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t, self.q)
        if self.stabilization is not None:
            self.stabilization.initializeElementQuadrature(self.mesh, self.timeIntegration.t, self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None:
            self.shockCapturing.initializeElementQuadrature(self.mesh, self.timeIntegration.t, self.q)

    def calculateElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on element boundaries.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.tmpvt.calculateElementBoundaryQuadrature()
        pass

    def calculateExteriorElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        if self.postProcessing:
            self.tmpvt.calculateExteriorElementBoundaryQuadrature()
        #
        # get physical locations of element boundary quadrature points
        #
        # assume all components live on the same mesh
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        self.fluxBoundaryConditionsObjectsDict = dict([(cj, FluxBoundaryConditions(self.mesh,
                                                                                   self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                   self.ebqe[('x')],
                                                                                   self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                   self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in list(self.advectiveFluxBoundaryConditionsSetterDict.keys())])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t, self.ebqe)

    def estimate_mt(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        self.h_dof_sge[:] = self.u[0].dof
        self.u_dof_sge[:] = self.u[1].dof
        self.v_dof_sge[:] = self.u[2].dof
        # if self.postProcessing:
        #     from proteus.cSW2D import calculateVelocityAverage as cva
        #     cva(self.mesh.nExteriorElementBoundaries_global,
        #         self.mesh.exteriorElementBoundariesArray,
        #         self.mesh.nInteriorElementBoundaries_global,
        #         self.mesh.interiorElementBoundariesArray,
        #         self.mesh.elementBoundaryElementsArray,
        #         self.mesh.elementBoundaryLocalElementBoundariesArray,
        #         self.u[1].femSpace.dofMap.l2g,
        #         self.u[1].dof,
        #         self.u[2].dof,
        #         self.u[3].dof,
        #         self.ebq[('v',0)],
        #         self.ebqe[('velocity',0)],
        #         self.ebq_global[('velocityAverage',0)])
        # self.sw2d.calculateVelocityAverage(self.mesh.nExteriorElementBoundaries_global,
        #                                      self.mesh.exteriorElementBoundariesArray,
        #                                      self.mesh.nInteriorElementBoundaries_global,
        #                                      self.mesh.interiorElementBoundariesArray,
        #                                      self.mesh.elementBoundaryElementsArray,
        #                                      self.mesh.elementBoundaryLocalElementBoundariesArray,
        #                                      self.mesh.nodeArray,
        #                                      self.mesh.elementNodesArray,
        #                                      self.u[0].femSpace.elementMaps.psi_trace,
        #                                      self.u[0].femSpace.elementMaps.grad_psi_trace,
        #                                      self.u[0].femSpace.elementMaps.boundaryNormals,
        #                                      self.u[0].femSpace.elementMaps.boundaryJacobians,
        #                                      self.u[1].femSpace.dofMap.l2g,
        #                                      self.u[1].dof,
        #                                      self.u[2].dof,
        #                                      self.u[3].dof,
        #                                      self.u[1].femSpace.psi_trace,
        #                                      self.ebqe[('velocity',0)],
        #                                      self.ebq_global[('velocityAverage',0)])
        OneLevelTransport.calculateAuxiliaryQuantitiesAfterStep(self)

    def getForce(self, cg, forceExtractionFaces, force, moment):
        pass
        # """
        # Calculate the element residuals and add in to the global residual
        # """

        # #Load the unknowns into the finite element dof
        # #self.timeIntegration.calculateCoefs()
        # #self.timeIntegration.calculateU(u)
        # #self.setUnknowns(self.timeIntegration.u)
        # #cek todo put in logic to skip if BC's don't depend on t or u
        # #hack
        # if self.bcsTimeDependent or not self.bcsSet:
        #     self.bcsSet=True
        #     #Dirichlet boundary conditions
        #     self.numericalFlux.setDirichletValues(self.ebqe)
        #     #Flux boundary conditions
        #     for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
        #         for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
        #             if self.coefficients.advection.has_key(ci):
        #                 self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
        #                 self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
        #         for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
        #             for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
        #                 self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
        #                 self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1

        # # Tag boundaries for force/moment extraction
        # #for cj in range(len(self.dirichletConditionsForceDOF)):
        # #     for dofN,g in self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.iteritems():

        # print "SW2D Force Faces",len(forceExtractionFaces)

        # #force  = numpy.zeros(3,'d')
        # #moment = numpy.zeros(3,'d')

        # self.Ct_sge = 4.0
        # self.Cd_sge = 144.0
        # self.C_b    = 10.0

        # self.sw2d.calculateForce(#element
        #     self.u[0].femSpace.elementMaps.psi,
        #     self.u[0].femSpace.elementMaps.grad_psi,
        #     self.mesh.nodeArray,
        #     self.mesh.elementNodesArray,
        #     self.elementQuadratureWeights[('u',0)],
        #     self.u[0].femSpace.psi,
        #     self.u[0].femSpace.grad_psi,
        #     self.u[0].femSpace.psi,
        #     self.u[0].femSpace.grad_psi,
        #     self.u[1].femSpace.psi,
        #     self.u[1].femSpace.grad_psi,
        #     self.u[1].femSpace.psi,
        #     self.u[1].femSpace.grad_psi,
        #     #element boundary
        #     self.u[0].femSpace.elementMaps.psi_trace,
        #     self.u[0].femSpace.elementMaps.grad_psi_trace,
        #     self.elementBoundaryQuadratureWeights[('u',0)],
        #     self.u[0].femSpace.psi_trace,
        #     self.u[0].femSpace.grad_psi_trace,
        #     self.u[0].femSpace.psi_trace,
        #     self.u[0].femSpace.grad_psi_trace,
        #     self.u[1].femSpace.psi_trace,
        #     self.u[1].femSpace.grad_psi_trace,
        #     self.u[1].femSpace.psi_trace,
        #     self.u[1].femSpace.grad_psi_trace,
        #     self.u[0].femSpace.elementMaps.boundaryNormals,
        #     self.u[0].femSpace.elementMaps.boundaryJacobians,
        #     #physics
        #     self.mesh.elementDiametersArray,
        #     self.stabilization.hFactor,
        #     self.mesh.nElements_global,
        #     self.coefficients.useRBLES,
        #     self.coefficients.useMetrics,
        #     self.timeIntegration.alpha_bdf,
        #     self.coefficients.epsFact_density,
        #     self.coefficients.epsFact,
        #     self.coefficients.sigma,
        #     self.coefficients.rho_0,
        #     self.coefficients.nu_0,
        #     self.coefficients.rho_1,
        #     self.coefficients.nu_1,
        #     self.Ct_sge,
        #     self.Cd_sge,
        #     self.shockCapturing.shockCapturingFactor,
        #     self.C_b,
        #     self.u[0].femSpace.dofMap.l2g,
        #     self.u[1].femSpace.dofMap.l2g,
        #     self.u[0].dof,
        #     self.u[1].dof,
        #     self.u[2].dof,
        #     self.u[3].dof,
        #     self.coefficients.g,
        #      self.q[('cfl',0)],   # ULTRA UGLY HACK self.q[('rho_0')],
        #     self.coefficients.q_phi,
        #     self.coefficients.q_n,
        #     self.coefficients.q_kappa,
        #     self.timeIntegration.m_tmp[1],
        #     self.timeIntegration.m_tmp[2],
        #     self.timeIntegration.m_tmp[3],
        #     self.q[('f',0)],
        #     self.timeIntegration.beta_bdf[1],
        #     self.timeIntegration.beta_bdf[2],
        #     self.timeIntegration.beta_bdf[3],
        #     self.stabilization.v_last,
        #     self.q[('cfl',0)],
        #     self.q[('numDiff',1,1)],
        #     self.q[('numDiff',2,2)],
        #     self.q[('numDiff',3,3)],
        #     self.shockCapturing.numDiff_last[1],
        #     self.shockCapturing.numDiff_last[2],
        #     self.shockCapturing.numDiff_last[3],
        #     self.coefficients.sdInfo[(1,1)][0],self.coefficients.sdInfo[(1,1)][1],
        #     self.coefficients.sdInfo[(1,2)][0],self.coefficients.sdInfo[(1,2)][1],
        #     self.coefficients.sdInfo[(1,3)][0],self.coefficients.sdInfo[(1,3)][1],
        #     self.coefficients.sdInfo[(2,2)][0],self.coefficients.sdInfo[(2,2)][1],
        #     self.coefficients.sdInfo[(2,1)][0],self.coefficients.sdInfo[(2,1)][1],
        #     self.coefficients.sdInfo[(2,3)][0],self.coefficients.sdInfo[(2,3)][1],
        #     self.coefficients.sdInfo[(3,3)][0],self.coefficients.sdInfo[(3,3)][1],
        #     self.coefficients.sdInfo[(3,1)][0],self.coefficients.sdInfo[(3,1)][1],
        #     self.coefficients.sdInfo[(3,2)][0],self.coefficients.sdInfo[(3,2)][1],
        #     self.offset[0],self.offset[1],self.offset[2],self.offset[3],
        #     self.stride[0],self.stride[1],self.stride[2],self.stride[3],
        #     cg, force, moment,
        #     self.mesh.nExteriorElementBoundaries_global,
        #     self.mesh.exteriorElementBoundariesArray,
        #     self.mesh.elementBoundaryElementsArray,
        #     self.mesh.elementBoundaryLocalElementBoundariesArray,
        #     forceExtractionFaces,len(forceExtractionFaces),
        #     self.coefficients.ebqe_phi,
        #     self.coefficients.ebqe_n,
        #     self.coefficients.ebqe_kappa,
        #     self.numericalFlux.isDOFBoundary[0],
        #     self.numericalFlux.isDOFBoundary[1],
        #     self.numericalFlux.isDOFBoundary[2],
        #     self.numericalFlux.isDOFBoundary[3],
        #     self.ebqe[('advectiveFlux_bc_flag',0)],
        #     self.ebqe[('advectiveFlux_bc_flag',1)],
        #     self.ebqe[('advectiveFlux_bc_flag',2)],
        #     self.ebqe[('advectiveFlux_bc_flag',3)],
        #     self.ebqe[('diffusiveFlux_bc_flag',1,1)],
        #     self.ebqe[('diffusiveFlux_bc_flag',2,2)],
        #     self.ebqe[('diffusiveFlux_bc_flag',3,3)],
        #     self.numericalFlux.ebqe[('u',0)],
        #     self.ebqe[('advectiveFlux_bc',0)],
        #     self.ebqe[('advectiveFlux_bc',1)],
        #     self.ebqe[('advectiveFlux_bc',2)],
        #     self.ebqe[('advectiveFlux_bc',3)],
        #     self.numericalFlux.ebqe[('u',1)],
        #     self.ebqe[('diffusiveFlux_bc',1,1)],
        #     self.ebqe[('penalty')],
        #     self.numericalFlux.ebqe[('u',2)],
        #     self.ebqe[('diffusiveFlux_bc',2,2)],
        #     self.numericalFlux.ebqe[('u',3)],
        #     self.ebqe[('diffusiveFlux_bc',3,3)],
        #     self.q[('velocity',0)],
        #     self.ebqe[('velocity',0)],
        #     self.ebq_global[('totalFlux',0)],
        #     self.elementResidual[0])

        # #from mpi4py import MPI
        # #comm = MPI.COMM_WORLD

        # #tmp1 = numpy.zeros(3,'d')
        # #tmp2 = numpy.zeros(3,'d')
        # #comm.Allreduce(force,  tmp1, op=MPI.SUM)
        # #comm.Allreduce(moment, tmp2, op=MPI.SUM)
        # #force  [:] = tmp1
        # #moment [:] = tmp2

        # for i in range(3):
        #       force[i]  = globalSum(force[i])
        #       moment[i] = globalSum(moment[i])

        # #simport time
        # #time.sleep(1)
        # ##comm.Barrier()
        # ##if self.comm.rank() == 0:
        # #print cg
        # #print "Force and moment in sw2d getForce"
        # #print force
        # #print moment
        # ##comm.Barrier()
        # #import time
        # #time.sleep(1)
