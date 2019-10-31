from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import proteus
from proteus.mprans.cGN_SW2DCV import *


class NumericalFlux(proteus.NumericalFlux.ShallowWater_2D):
    hasInterior = False

    def __init__(self, vt, getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None,
                 h_eps=1.0e-5,
                 tol_u=1.0e-5):
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


class RKEV(proteus.TimeIntegration.SSP):
    from proteus import TimeIntegration
    """
    Wrapper for SSPRK time integration using EV

    ... more to come ...
    """

    def __init__(self, transport, timeOrder=1, runCFL=0.1, integrateInterpolationPoints=False):
        BackwardEuler.__init__(
            self, transport, integrateInterpolationPoints=integrateInterpolationPoints)
        self.trasport = transport
        self.runCFL = runCFL
        self.dtLast = None
        self.dtRatioMax = 2.0
        self.isAdaptive = True
        # About the cfl
        assert hasattr(
            transport, 'edge_based_cfl'), "No edge based cfl defined"
        self.edge_based_cfl = transport.edge_based_cfl
        # Stuff particular for SSP
        self.timeOrder = timeOrder  # order of approximation
        self.nStages = timeOrder  # number of stages total
        self.lstage = 0  # last stage completed
        # storage vectors (at old time step)
        self.u_dof_last = {}
        # per component lstage values
        self.u_dof_lstage = {}
        for ci in range(self.nc):
            self.u_dof_last[ci] = transport.u[ci].dof.copy()
            self.u_dof_lstage[ci] = transport.u[ci].dof.copy()

    def choose_dt(self):
        maxCFL = 1.0e-6
        # COMPUTE edge_based_cfl
        rowptr_cMatrix, colind_cMatrix, Cx = self.transport.cterm_global[0].getCSRrepresentation(
        )
        rowptr_cMatrix, colind_cMatrix, Cy = self.transport.cterm_global[1].getCSRrepresentation(
        )
        rowptr_cMatrix, colind_cMatrix, CTx = self.transport.cterm_global_transpose[0].getCSRrepresentation(
        )
        rowptr_cMatrix, colind_cMatrix, CTy = self.transport.cterm_global_transpose[1].getCSRrepresentation(
        )
        numDOFsPerEqn = self.transport.u[0].dof.size

        adjusted_maxCFL = self.transport.sw2d.calculateEdgeBasedCFL(
            self.transport.coefficients.g,
            numDOFsPerEqn,
            self.transport.ML,
            self.transport.u[0].dof,
            self.transport.u[1].dof,
            self.transport.u[2].dof,
            self.transport.u[3].dof,
            self.transport.coefficients.b.dof,
            rowptr_cMatrix,
            colind_cMatrix,
            self.transport.hEps,
            self.transport.hReg,
            Cx,
            Cy,
            CTx,
            CTy,
            self.transport.dLow,
            self.runCFL,
            self.transport.edge_based_cfl,
            0)

        maxCFL = max(maxCFL, max(adjusted_maxCFL,
                                 globalMax(self.edge_based_cfl.max())))
        self.dt = old_div(self.runCFL, maxCFL)

        if self.dtLast is None:
            self.dtLast = self.dt
        if old_div(self.dt, self.dtLast) > self.dtRatioMax:
            self.dt = self.dtLast * self.dtRatioMax
        self.t = self.tLast + self.dt
        # Ignoring dif. time step levels
        self.substeps = [self.t for i in range(self.nStages)]

        assert (self.dt > 1E-8), ("Time step is probably getting too small: ", self.dt, adjusted_maxCFL,)


    def initialize_dt(self, t0, tOut, q):
        """
        Modify self.dt
        """
        self.tLast = t0
        self.choose_dt()
        self.t = t0 + self.dt

    def setCoefficients(self):
        """
        beta are all 1's here
        mwf not used right now
        """
        # Not needed for an implementation when alpha and beta are not used

    def updateStage(self):
        """
        Need to switch to use coefficients
        """
        self.lstage += 1
        assert self.timeOrder in [1, 2, 3]
        assert self.lstage > 0 and self.lstage <= self.timeOrder
        # print "within update stage...: ", self.lstage
        if self.timeOrder == 3:
            if self.lstage == 1:
                logEvent("First stage of SSP33 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                # update u_dof_old
                self.transport.h_dof_old[:] = self.u_dof_lstage[0]
                self.transport.hu_dof_old[:] = self.u_dof_lstage[1]
                self.transport.hv_dof_old[:] = self.u_dof_lstage[2]
                self.transport.heta_dof_old[:] = self.u_dof_lstage[3]
                self.transport.hw_dof_old[:] = self.u_dof_lstage[4]
            elif self.lstage == 2:
                logEvent("Second stage of SSP33 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                    self.u_dof_lstage[ci] *= old_div(1., 4.)
                    self.u_dof_lstage[ci] += 3. / 4. * self.u_dof_last[ci]
                # update u_dof_old
                self.transport.h_dof_old[:] = self.u_dof_lstage[0]
                self.transport.hu_dof_old[:] = self.u_dof_lstage[1]
                self.transport.hv_dof_old[:] = self.u_dof_lstage[2]
                self.transport.heta_dof_old[:] = self.u_dof_lstage[3]
                self.transport.hw_dof_old[:] = self.u_dof_lstage[4]
            else:
                logEvent("Third stage of SSP33 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                    self.u_dof_lstage[ci][:] *= old_div(2.0, 3.0)
                    self.u_dof_lstage[ci][:] += 1.0 / 3.0 * self.u_dof_last[ci]
                    # update solution to u[0].dof
                    self.transport.u[ci].dof[:] = self.u_dof_lstage[ci]
                # update u_dof_old
                self.transport.h_dof_old[:] = self.u_dof_last[0]
                self.transport.hu_dof_old[:] = self.u_dof_last[1]
                self.transport.hv_dof_old[:] = self.u_dof_last[2]
                self.transport.heta_dof_old[:] = self.u_dof_last[3]
                self.transport.hw_dof_old[:] = self.u_dof_last[4]
        elif self.timeOrder == 2:
            if self.lstage == 1:
                logEvent("First stage of SSP22 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                # Update u_dof_old
                self.transport.h_dof_old[:] = self.u_dof_lstage[0]
                self.transport.hu_dof_old[:] = self.u_dof_lstage[1]
                self.transport.hv_dof_old[:] = self.u_dof_lstage[2]
                self.transport.heta_dof_old[:] = self.u_dof_lstage[3]
                self.transport.hw_dof_old[:] = self.u_dof_lstage[4]
            else:
                logEvent("Second stage of SSP22 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                    self.u_dof_lstage[ci][:] *= old_div(1., 2.)
                    self.u_dof_lstage[ci][:] += 1. / 2. * self.u_dof_last[ci]
                    # update solution to u[0].dof
                    self.transport.u[ci].dof[:] = self.u_dof_lstage[ci]
                # Update u_dof_old
                self.transport.h_dof_old[:] = self.u_dof_last[0]
                self.transport.hu_dof_old[:] = self.u_dof_last[1]
                self.transport.hv_dof_old[:] = self.u_dof_last[2]
                self.transport.heta_dof_old[:] = self.u_dof_last[3]
                self.transport.hw_dof_old[:] = self.u_dof_last[4]
        else:
            assert self.timeOrder == 1
            logEvent("FE method finished", level=4)

    def initializeTimeHistory(self, resetFromDOF=True):
        """
        Push necessary information into time history arrays
        """
        for ci in range(self.nc):
            self.u_dof_last[ci][:] = self.transport.u[ci].dof[:]

    def updateTimeHistory(self, resetFromDOF=False):
        """
        assumes successful step has been taken
        """
        self.t = self.tLast + self.dt
        for ci in range(self.nc):
            self.u_dof_last[ci][:] = self.transport.u[ci].dof[:]
        self.lstage = 0
        self.dtLast = self.dt
        self.tLast = self.t

    def generateSubsteps(self, tList):
        """
        create list of substeps over time values given in tList. These correspond to stages
        """
        self.substeps = []
        tLast = self.tLast
        for t in tList:
            dttmp = t - tLast
            self.substeps.extend([tLast + dttmp for i in range(self.nStages)])
            tLast = t

    def resetOrder(self, order):
        """
        initialize data structures for stage updges
        """
        self.timeOrder = order  # order of approximation
        self.nStages = order  # number of stages total
        self.lstage = 0  # last stage completed
        # storage vectors
        # per component stage values, list with array at each stage
        self.u_dof_lstage = {}
        for ci in range(self.nc):
            self.u_dof_lstage[ci] = self.transport.u[ci].dof.copy()
        self.substeps = [self.t for i in range(self.nStages)]

    def setFromOptions(self, nOptions):
        """
        allow classes to set various numerical parameters
        """
        if 'runCFL' in dir(nOptions):
            self.runCFL = nOptions.runCFL
        flags = ['timeOrder']
        for flag in flags:
            if flag in dir(nOptions):
                val = getattr(nOptions, flag)
                setattr(self, flag, val)
                if flag == 'timeOrder':
                    self.resetOrder(self.timeOrder)


class Coefficients(proteus.TransportCoefficients.TC_base):
    """
    The coefficients for the shallow water, modified Green-Naghdi equations
    """

    def __init__(self,
                 bathymetry,
                 nu=1.004e-6,
                 g=9.8,
                 nd=2,
                 sd=True,
                 movingDomain=False,
                 useRBLES=0.0,
                 useMetrics=0.0,
                 modelIndex=0,
                 cE=1.0,
                 LUMPED_MASS_MATRIX=1,
                 LINEAR_FRICTION=0,
                 mannings=0.,
                 forceStrongConditions=True,
                 constrainedDOFs=None):
        self.forceStrongConditions = forceStrongConditions
        self.constrainedDOFs = constrainedDOFs
        self.bathymetry = bathymetry
        self.useRBLES = useRBLES
        self.useMetrics = useMetrics
        self.sd = sd
        self.nu = nu
        self.g = g
        self.nd = nd
        self.cE = cE
        self.LUMPED_MASS_MATRIX = LUMPED_MASS_MATRIX
        self.LINEAR_FRICTION = LINEAR_FRICTION
        self.mannings = mannings
        self.modelIndex = modelIndex
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction = {}
        hamiltonian = {}
        if nd == 2:
            variableNames = ['h', 'h_u', 'h_v', 'h_eta', 'h_w']
            mass = {0: {0: 'linear'},
                    1: {0: 'linear', 1: 'linear'},
                    2: {0: 'linear', 2: 'linear'},
                    3: {0: 'linear', 3: 'linear'},
                    4: {0: 'linear', 4: 'linear'}}
            advection = {0: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear',
                             4: 'nonlinear'},
                         1: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear',
                             4: 'nonlinear'},
                         2: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear',
                             4: 'nonlinear'},
                         3: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear',
                             4: 'nonlinear'},
                         4: {0: 'nonlinear',
                             1: 'nonlinear',
                             2: 'nonlinear',
                             3: 'nonlinear',
                             4: 'nonlinear'}}
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
                             5,  # Number of components
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
        self.model = modelList[self.modelIndex]
        # pass

    def initializeMesh(self, mesh):
        x = mesh.nodeArray[:, 0]
        y = mesh.nodeArray[:, 1]
        if self.bathymetry is None:
            self.b.dof = mesh.nodeArray[:, 2].copy()
        else:
            self.b.dof = self.bathymetry[0]([x, y])
        # self.b.dof[:] = 0. #TMP

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

    def preStep(self, t, firstStep=False):
        if firstStep:
            # Init boundaryIndex
            assert self.model.boundaryIndex is None and self.model.normalx is not None , \
                    "Check boundaryIndex, normalx and normaly"
            self.model.boundaryIndex = []
            for i in range(self.model.normalx.size):
                if self.model.normalx[i] != 0 or self.model.normaly[i] != 0:
                    self.model.boundaryIndex.append(i)
            self.model.boundaryIndex = np.array(self.model.boundaryIndex)
        #
        self.model.h_dof_old[:] = self.model.u[0].dof
        self.model.hu_dof_old[:] = self.model.u[1].dof
        self.model.hv_dof_old[:] = self.model.u[2].dof
        self.model.heta_dof_old[:] = self.model.u[3].dof
        self.model.hw_dof_old[:] = self.model.u[4].dof

    def postStep(self, t, firstStep=False):
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
                 name='GN_SW2DCV',
                 reuse_trial_and_test_quadrature=True,
                 sd=True,
                 movingDomain=False,
                 bdyNullSpace=False):
        self.bdyNullSpace = bdyNullSpace
        self.inf_norm_hu = []  # To test 1D well balancing
        self.secondCallCalculateResidual = 0
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
                assert self.u[ci].femSpace.__class__.__name__ == self.u[
                    0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        # Simplicial Mesh
        # assume the same mesh for  all components for now
        self.mesh = self.u[0].femSpace.mesh
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        # explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.dirichletNodeSetList = None
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
        # no velocity post-processing for now
        self.conservativeFlux = conservativeFluxDict
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        # determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci in range(self.nc):
            self.elementBoundaryIntegrals[ci] = ((self.conservativeFlux is not None)
                                                 or (numericalFluxType is not None)
                                                 or (self.fluxBoundaryConditions[ci] == 'outFlow')
                                                 or (self.fluxBoundaryConditions[ci] == 'mixedFlow')
                                                 or (self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        # calculate some dimensions
        #
        # assume same space dim for all variables
        self.nSpace_global = self.u[0].femSpace.nSpace_global
        self.nDOF_trial_element = [
            u_j.femSpace.max_nDOF_element for u_j in list(self.u.values())]
        self.nDOF_phi_trial_element = [
            phi_k.femSpace.max_nDOF_element for phi_k in list(self.phi.values())]
        self.n_phi_ip_element = [
            phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for phi_k in list(self.phi.values())]
        self.nDOF_test_element = [
            femSpace.max_nDOF_element for femSpace in list(self.testSpace.values())]
        self.nFreeDOF_global = [dc.nFreeDOF_global for dc in list(
            self.dirichletConditions.values())]
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
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if ('numDiff', ci, ci) in elementQuadrature:
                        elementQuadratureDict[(
                            'numDiff', ci, ci)] = elementQuadrature[('numDiff', ci, ci)]
                    else:
                        elementQuadratureDict[(
                            'numDiff', ci, ci)] = elementQuadrature['default']
                else:
                    elementQuadratureDict[(
                        'numDiff', ci, ci)] = elementQuadrature
        if massLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('m', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[(
                    'stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        if reactionLumping:
            for ci in list(self.coefficients.mass.keys()):
                elementQuadratureDict[('r', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[(
                    'stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
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
        self.nQuadraturePoints_global = self.nQuadraturePoints_element * \
            self.mesh.nElements_global
        #
        # Repeat the same thing for the element boundary quadrature
        #
        (self.elementBoundaryQuadraturePoints,
         self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[
            0]
        self.nElementBoundaryQuadraturePoints_global = (self.mesh.nElements_global
                                                        * self.mesh.nElementBoundaries_element
                                                        * self.nElementBoundaryQuadraturePoints_elementBoundary)

        #
        # simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q = {}
        self.ebq = {}
        self.ebq_global = {}
        self.ebqe = {}
        self.phi_ip = {}
        # To compute edge_based_cfl from withing choose_dt of RKEV
        self.edge_based_cfl = numpy.zeros(self.u[0].dof.shape)
        self.dLow = None
        self.hBT = None
        self.huBT = None
        self.hvBT = None
        self.hetaBT = None
        self.hwBT = None
        # Old DOFs
        # NOTE (Mql): It is important to link h_dof_old by reference with u[0].dof (and so on).
        # This is because  I need the initial condition to be passed to them as well (before calling calculateResidual).
        # During preStep I change this and copy the values instead of keeping the reference.
        self.h_dof_old = None
        self.hu_dof_old = None
        self.hv_dof_old = None
        self.heta_dof_old = None
        self.hw_dof_old = None

        # Vector for mass matrix
        self.check_positivity_water_height = True
        # mesh
        self.h_dof_sge = self.u[0].dof.copy()
        self.hu_dof_sge = self.u[1].dof.copy()
        self.hv_dof_sge = self.u[2].dof.copy()
        self.heta_dof_sge = self.u[3].dof.copy()
        self.hw_dof_sge = self.u[4].dof.copy()
        self.q['x'] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element, 3), 'd')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                      self.nElementBoundaryQuadraturePoints_elementBoundary, 3), 'd')
        self.ebq_global[('totalFlux', 0)] = numpy.zeros(
            (self.mesh.nElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebq_global[('velocityAverage', 0)] = numpy.zeros((self.mesh.nElementBoundaries_global,
                                                               self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.q[('dV_u', 0)] = (old_div(1.0, self.mesh.nElements_global)) * \
            numpy.ones((self.mesh.nElements_global,
                        self.nQuadraturePoints_element), 'd')
        self.q[('dV_u', 1)] = (old_div(1.0, self.mesh.nElements_global)) * \
            numpy.ones((self.mesh.nElements_global,
                        self.nQuadraturePoints_element), 'd')
        self.q[('dV_u', 2)] = (old_div(1.0, self.mesh.nElements_global)) * \
            numpy.ones((self.mesh.nElements_global,
                        self.nQuadraturePoints_element), 'd')
        self.q[('u', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 3)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('u', 4)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m', 0)] = self.q[('u', 0)]
        self.q[('m', 1)] = self.q[('u', 1)]
        self.q[('m', 2)] = self.q[('u', 2)]
        self.q[('m', 3)] = self.q[('u', 3)]
        self.q[('m', 4)] = self.q[('u', 4)]
        self.q[('m_last', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_last', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 1)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 2)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('f', 0)] = numpy.zeros((self.mesh.nElements_global,
                                        self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('velocity', 0)] = numpy.zeros((self.mesh.nElements_global,
                                               self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('cfl', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.ebqe[('u', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 3)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('u', 4)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                           self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc_flag', 0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 3)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc_flag', 4)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 1, 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 2, 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 3, 3)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('diffusiveFlux_bc_flag', 4, 4)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc', 0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 3)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux_bc', 4)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 1, 1)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('penalty')] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                              self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 2, 2)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 3, 3)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('diffusiveFlux_bc', 4, 4)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('velocity', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('velocity', 1)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('velocity', 2)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                  self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set(
            [('u', ci) for ci in range(self.nc)])
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
        logEvent(
            "Interpolation points for nonlinear diffusion potential (phi_ip)", level=9)
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
            self.inflowBoundaryBC[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,), 'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global, self.nDOF_trial_element[cj]), 'd')
            self.inflowFlux[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
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
        logEvent(memory("inflowBC, internalNodes,updateLocal2Global",
                        "OneLevelTransport"), level=4)
        self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        logEvent(memory("TimeIntegration", "OneLevelTransport"), level=4)
        logEvent("Calculating numerical quadrature formulas", 2)
        self.calculateQuadrature()
        self.setupFieldStrides()

        # hEps: this is use to regularize the flux and re-define the dry states
        self.hEps = None
        self.hReg = None
        self.ML = None  # lumped mass matrix
        self.MC_global = None  # consistent mass matrix
        # Global C Matrices (mql)
        self.cterm_global = None
        self.cterm_transpose_global = None
        # For FCT
        self.extendedSourceTerm_hu = None
        self.extendedSourceTerm_hv = None
        self.extendedSourceTerm_heta = None
        self.extendedSourceTerm_hw = None
        self.new_SourceTerm_hu = None
        self.new_SourceTerm_hv = None
        self.new_SourceTerm_heta = None
        self.new_SourceTerm_hw = None
        self.dH_minus_dL = None
        self.muH_minus_muL = None
        # NORMALS
        self.COMPUTE_NORMALS = 1
        self.normalx = None
        self.normaly = None
        self.boundaryIndex = None
        self.reflectingBoundaryConditions = False

        if 'reflecting_BCs' in dir(options) and options.reflecting_BCs == 1:
            self.reflectingBoundaryConditions = True
        # Aux quantity at DOFs to be filled by optimized code (MQL)
        self.quantDOFs = None

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
                    self.ebq_global['penalty'][ebN, k] = old_div(self.numericalFlux.penalty_constant,
                                                                 (self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power))
        # penalty term
        # cek move  to Numerical flux initialization
        if 'penalty' in self.ebqe:
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE, k] = old_div(self.numericalFlux.penalty_constant,
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
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(
            self)
        logEvent(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        for ci, fbcObject in list(self.fluxBoundaryConditionsObjectsDict.items()):
            self.ebqe[('advectiveFlux_bc_flag', ci)] = numpy.zeros(
                self.ebqe[('advectiveFlux_bc', ci)].shape, 'i')
            for t, g in list(fbcObject.advectiveFluxBoundaryConditionsDict.items()):
                if ci in self.coefficients.advection:
                    self.ebqe[('advectiveFlux_bc', ci)][t[0], t[1]] = g(
                        self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1
            for ck, diffusiveFluxBoundaryConditionsDict in list(fbcObject.diffusiveFluxBoundaryConditionsDictDict.items()):
                self.ebqe[('diffusiveFlux_bc_flag', ck, ci)] = numpy.zeros(
                    self.ebqe[('diffusiveFlux_bc', ck, ci)].shape, 'i')
                for t, g in list(diffusiveFluxBoundaryConditionsDict.items()):
                    self.ebqe[('diffusiveFlux_bc', ck, ci)][t[0], t[1]] = g(
                        self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag', ck, ci)
                              ][t[0], t[1]] = 1
        # self.numericalFlux.setDirichletValues(self.ebqe)
        if self.movingDomain:
            self.MOVING_DOMAIN = 1.0
        else:
            self.MOVING_DOMAIN = 0.0
        # cek hack
        self.movingDomain = False
        self.MOVING_DOMAIN = 0.0
        if self.mesh.nodeVelocityArray is None:
            self.mesh.nodeVelocityArray = numpy.zeros(
                self.mesh.nodeArray.shape, 'd')
        # cek/ido todo replace python loops in modules with optimized code if possible/necessary
        self.forceStrongConditions = self.coefficients.forceStrongConditions
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
        self.sw2d = cGN_SW2DCV_base(self.nSpace_global,
                                    self.nQuadraturePoints_element,
                                    self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                    self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                    self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                    self.nElementBoundaryQuadraturePoints_elementBoundary,
                                    compKernelFlag)

        self.calculateResidual = self.sw2d.calculateResidual
        if (self.coefficients.LUMPED_MASS_MATRIX):
            self.calculateJacobian = self.sw2d.calculateLumpedMassMatrix
            print("calculating lumped ")
        else:
            self.calculateJacobian = self.sw2d.calculateMassMatrix
            print("calculating full mass matrix")
        #
        self.dofsXCoord = None
        self.dofsYCoord = None
        self.constrainedDOFsIndices = None
        self.dataStructuresInitialized = False

    def FCTStep(self):
        # NOTE: this function is meant to be called within the solver
        rowptr, colind, MassMatrix = self.MC_global.getCSRrepresentation()
        # Extract hnp1 from global solution u
        index = list(range(0, len(self.timeIntegration.u)))
        hIndex = index[0::5]
        huIndex = index[1::5]
        hvIndex = index[2::5]
        hetaIndex = index[3::5]
        hwIndex = index[4::5]
        # create limited solution
    #     limited_hnp1 = numpy.zeros(self.h_dof_old.shape)
    #     limited_hunp1 = numpy.zeros(self.h_dof_old.shape)
    #     limited_hvnp1 = numpy.zeros(self.h_dof_old.shape)
    #     limited_hetanp1 = numpy.zeros(self.h_dof_old.shape)
    #     limited_hwnp1 = numpy.zeros(self.h_dof_old.shape)
    # #     # Do some type of limitation
    #     self.sw2d.convexLimiting(self.timeIntegration.dt,
    #                              # self.sw2d.FCTStep(self.timeIntegration.dt,
    #                              self.nnz,  # number of non zero entries
    #                              len(rowptr) - 1,  # number of DOFs
    #                              self.ML,  # Lumped mass matrix
    #                              self.h_dof_old,
    #                              self.hu_dof_old,
    #                              self.hv_dof_old,
    #                              self.heta_dof_old,
    #                              self.hw_dof_old,
    #                              self.coefficients.b.dof,
    #                              # high order solution
    #                              self.timeIntegration.u[hIndex],
    #                              self.timeIntegration.u[huIndex],
    #                              self.timeIntegration.u[hvIndex],
    #                              self.timeIntegration.u[hetaIndex],
    #                              self.timeIntegration.u[hwIndex],
    #                              self.extendedSourceTerm_hu,
    #                              self.extendedSourceTerm_hv,
    #                              self.extendedSourceTerm_heta,
    #                              self.extendedSourceTerm_hw,
    #                              limited_hnp1,
    #                              limited_hunp1,
    #                              limited_hvnp1,
    #                              limited_hetanp1,
    #                              limited_hwnp1,
    #                              # Row indices for Sparsity Pattern (convenient for DOF loops)
    #                              rowptr,
    #                              # Column indices for Sparsity Pattern (convenient for DOF loops)
    #                              colind,
    #                              MassMatrix,
    #                              self.dH_minus_dL,
    #                              self.muH_minus_muL,
    #                              self.hEps,
    #                              self.hReg,
    #                              self.coefficients.LUMPED_MASS_MATRIX,
    #                              self.dLow,
    #                              self.hBT,
    #                              self.huBT,
    #                              self.hvBT,
    #                              self.hetaBT,
    #                              self.hwBT,
    #                              self.new_SourceTerm_hu,
    #                              self.new_SourceTerm_hv,
    #                              self.new_SourceTerm_heta,
    #                              self.new_SourceTerm_hw)
    #
    #     # Pass the post processed hnp1 solution to global solution u
    #     self.timeIntegration.u[hIndex] = limited_hnp1
    #     self.timeIntegration.u[huIndex] = limited_hunp1
    #     self.timeIntegration.u[hvIndex] = limited_hvnp1
    #     self.timeIntegration.u[hetaIndex] = limited_hetanp1
    #     self.timeIntegration.u[hwIndex] = limited_hwnp1

    def getDOFsCoord(self):
        # get x,y coordinates of all DOFs #
        self.dofsXCoord = numpy.zeros(self.u[0].dof.shape, 'd')
        self.dofsYCoord = numpy.zeros(self.u[0].dof.shape, 'd')
        self.dofsXCoord[:] = self.mesh.nodeArray[:, 0]
        self.dofsYCoord[:] = self.mesh.nodeArray[:, 1]
        #

    def getCMatrices(self):
        # since we only need cterm_global to persist, we can drop the other self.'s
        self.cterm = {}
        self.cterm_a = {}
        self.cterm_global = {}
        self.cterm_transpose = {}
        self.cterm_a_transpose = {}
        self.cterm_global_transpose = {}
        # Sparsity pattern for Jacobian
        rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
        nnz = nzval.shape[-1]  # number of non-zero entries in sparse matrix

        ###########################################
        ##### SPARSITY PATTERN FOR C MATRICES #####
        ###########################################
        # construct nnz_cMatrix, czval_cMatrix, rowptr_cMatrix, colind_cMatrix C matrix
        nnz_cMatrix = nnz // 5 // 5  # This is always true for the modified Green-Naghdi in 2D
        nzval_cMatrix = numpy.zeros(nnz_cMatrix)
        rowptr_cMatrix = numpy.zeros(self.u[0].dof.size + 1, 'i')
        colind_cMatrix = numpy.zeros(nnz_cMatrix, 'i')
        # fill vector rowptr_cMatrix
        for i in range(1, rowptr_cMatrix.size):
            rowptr_cMatrix[i] = rowptr_cMatrix[i - 1] + \
                old_div((rowptr[5 * (i - 1) + 1] - rowptr[5 * (i - 1)]), 5)
            # = rowptr_cMatrix[i-1] + 1/3*(Number of columns of Jacobian's row 3*(i-1)=0, 3, 6, 9, 12, ... 3*(i-1), ..., 3*n-3)

        # fill vector colind_cMatrix
        i_cMatrix = 0  # ith row of cMatrix
        # 0 to total num of DOFs (i.e. num of rows of jacobian)
        for i in range(rowptr.size - 1):
            if (i % 5 == 0):  # Just consider the rows related to the h variable
                for j, offset in enumerate(range(rowptr[i], rowptr[i + 1])):
                    offset_cMatrix = list(
                        range(rowptr_cMatrix[i_cMatrix], rowptr_cMatrix[i_cMatrix + 1]))
                    if (j % 5 == 0):
                        colind_cMatrix[offset_cMatrix[old_div(j, 5)]] = old_div(
                            colind[offset], 5)
                i_cMatrix += 1
        # END OF SPARSITY PATTERN FOR C MATRICES

        di = numpy.zeros((self.mesh.nElements_global,
                          self.nQuadraturePoints_element,
                          self.nSpace_global),
                         'd')  # direction of derivative
        # JACOBIANS (FOR ELEMENT TRANSFORMATION)
        self.q[('J')] = numpy.zeros((self.mesh.nElements_global,
                                     self.nQuadraturePoints_element,
                                     self.nSpace_global,
                                     self.nSpace_global),
                                    'd')
        self.q[('inverse(J)')] = numpy.zeros((self.mesh.nElements_global,
                                              self.nQuadraturePoints_element,
                                              self.nSpace_global,
                                              self.nSpace_global),
                                             'd')
        self.q[('det(J)')] = numpy.zeros((self.mesh.nElements_global,
                                          self.nQuadraturePoints_element),
                                         'd')
        self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                         self.q['J'],
                                                         self.q['inverse(J)'],
                                                         self.q['det(J)'])
        self.q['abs(det(J))'] = numpy.abs(self.q['det(J)'])
        # SHAPE FUNCTIONS
        self.q[('w', 0)] = numpy.zeros((self.mesh.nElements_global,
                                        self.nQuadraturePoints_element,
                                        self.nDOF_test_element[0]),
                                       'd')
        self.q[('w*dV_m', 0)] = self.q[('w', 0)].copy()
        self.u[0].femSpace.getBasisValues(
            self.elementQuadraturePoints, self.q[('w', 0)])
        cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u', 0)],
                                             self.q['abs(det(J))'],
                                             self.q[('w', 0)],
                                             self.q[('w*dV_m', 0)])
        # GRADIENT OF TEST FUNCTIONS
        self.q[('grad(w)', 0)] = numpy.zeros((self.mesh.nElements_global,
                                              self.nQuadraturePoints_element,
                                              self.nDOF_test_element[0],
                                              self.nSpace_global),
                                             'd')
        self.u[0].femSpace.getBasisGradientValues(self.elementQuadraturePoints,
                                                  self.q['inverse(J)'],
                                                  self.q[('grad(w)', 0)])
        self.q[('grad(w)*dV_f', 0)] = numpy.zeros((self.mesh.nElements_global,
                                                   self.nQuadraturePoints_element,
                                                   self.nDOF_test_element[0],
                                                   self.nSpace_global),
                                                  'd')
        cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('u', 0)],
                                                      self.q['abs(det(J))'],
                                                      self.q[('grad(w)', 0)],
                                                      self.q[('grad(w)*dV_f', 0)])
        #
        # lumped mass matrix
        #
        # assume a linear mass term
        dm = np.ones(self.q[('u', 0)].shape, 'd')
        elementMassMatrix = np.zeros((self.mesh.nElements_global,
                                      self.nDOF_test_element[0],
                                      self.nDOF_trial_element[0]), 'd')
        cfemIntegrals.updateMassJacobian_weak_lowmem(dm,
                                                     self.q[('w', 0)],
                                                     self.q[('w*dV_m', 0)],
                                                     elementMassMatrix)
        self.MC_a = nzval_cMatrix.copy()
        self.MC_global = SparseMat(self.nFreeDOF_global[0],
                                   self.nFreeDOF_global[0],
                                   nnz_cMatrix,
                                   self.MC_a,
                                   colind_cMatrix,
                                   rowptr_cMatrix)
        cfemIntegrals.zeroJacobian_CSR(nnz_cMatrix, self.MC_global)
        cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                  self.l2g[0]['freeLocal'],
                                                                  self.l2g[0]['nFreeDOF'],
                                                                  self.l2g[0]['freeLocal'],
                                                                  self.csrRowIndeces[(
                                                                      0, 0)] // 5 // 5,
                                                                  old_div(
                                                                      self.csrColumnOffsets[(0, 0)], 5),
                                                                  elementMassMatrix,
                                                                  self.MC_global)
        diamD2 = numpy.sum(self.q['abs(det(J))'][:]
                           * self.elementQuadratureWeights[('u', 0)])
        self.ML = np.zeros((self.nFreeDOF_global[0],), 'd')
        self.hReg = np.zeros((self.nFreeDOF_global[0],), 'd')
        for i in range(self.nFreeDOF_global[0]):
            self.ML[i] = self.MC_a[rowptr_cMatrix[i]:rowptr_cMatrix[i + 1]].sum()
            self.hReg[i] = self.ML[i] / diamD2 * self.u[0].dof.max()
        # np.testing.assert_almost_equal(self.ML.sum(), self.mesh.volume, err_msg="Trace of lumped mass matrix should be the domain volume",verbose=True)
        # np.testing.assert_almost_equal(self.ML.sum(), diamD2, err_msg="Trace of lumped mass matrix should be the domain volume",verbose=True)

        for d in range(self.nSpace_global):  # spatial dimensions
            # C matrices
            self.cterm[d] = numpy.zeros((self.mesh.nElements_global,
                                         self.nDOF_test_element[0],
                                         self.nDOF_trial_element[0]), 'd')
            self.cterm_a[d] = nzval_cMatrix.copy()
            self.cterm_global[d] = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                                self.nFreeDOF_global[0],
                                                                nnz_cMatrix,
                                                                self.cterm_a[d],
                                                                colind_cMatrix,
                                                                rowptr_cMatrix)
            cfemIntegrals.zeroJacobian_CSR(nnz_cMatrix, self.cterm_global[d])
            di[:] = 0.0
            di[..., d] = 1.0
            cfemIntegrals.updateHamiltonianJacobian_weak_lowmem(di,
                                                                self.q[(
                                                                    'grad(w)*dV_f', 0)],
                                                                self.q[(
                                                                    'w', 0)],
                                                                self.cterm[d])  # int[(di*grad(wj))*wi*dV]
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(
                                                                          0, 0)] // 5 // 5,
                                                                      old_div(
                                                                          self.csrColumnOffsets[(0, 0)], 5),
                                                                      self.cterm[d],
                                                                      self.cterm_global[d])
            # C Transpose matrices
            self.cterm_transpose[d] = numpy.zeros((self.mesh.nElements_global,
                                                   self.nDOF_test_element[0],
                                                   self.nDOF_trial_element[0]), 'd')
            self.cterm_a_transpose[d] = nzval_cMatrix.copy()
            self.cterm_global_transpose[d] = LinearAlgebraTools.SparseMat(self.nFreeDOF_global[0],
                                                                          self.nFreeDOF_global[0],
                                                                          nnz_cMatrix,
                                                                          self.cterm_a_transpose[d],
                                                                          colind_cMatrix,
                                                                          rowptr_cMatrix)
            cfemIntegrals.zeroJacobian_CSR(
                nnz_cMatrix, self.cterm_global_transpose[d])
            di[:] = 0.0
            di[..., d] = -1.0
            cfemIntegrals.updateAdvectionJacobian_weak_lowmem(di,
                                                              self.q[('w', 0)],
                                                              self.q[(
                                                                  'grad(w)*dV_f', 0)],
                                                              self.cterm_transpose[d])  # -int[(-di*grad(wi))*wj*dV]
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(
                                                                          0, 0)] // 5 // 5,
                                                                      old_div(
                                                                          self.csrColumnOffsets[(0, 0)], 5),
                                                                      self.cterm_transpose[d],
                                                                      self.cterm_global_transpose[d])
        #
        self.rowptr_cMatrix, self.colind_cMatrix, self.Cx = self.cterm_global[0].getCSRrepresentation(
        )
        rowptr_cMatrix, colind_cMatrix, self.Cy = self.cterm_global[1].getCSRrepresentation(
        )
        rowptr_cMatrix, colind_cMatrix, self.CTx = self.cterm_global_transpose[0].getCSRrepresentation(
        )
        rowptr_cMatrix, colind_cMatrix, self.CTy = self.cterm_global_transpose[1].getCSRrepresentation(
        )
        # (mql): I am assuming all variables live on the same FE space
        self.numDOFsPerEqn = self.u[0].dof.size
        self.numNonZeroEntries = len(self.Cx)
    #

    def updateConstrainedDOFs(self):
        # get indices for constrained DOFs
        if self.constrainedDOFsIndices is None:
            self.constrainedDOFsIndices = []
            self.constrainedDOFsIndices = self.coefficients.constrainedDOFs[0](self.dofsXCoord,
                                                                               self.dofsYCoord)
        for i in self.constrainedDOFsIndices:
            x = self.dofsXCoord[i]
            y = self.dofsYCoord[i]
            (h, hu, hv, heta, hw) = self.coefficients.constrainedDOFs[1](x, y, self.timeIntegration.t,
                                                                         self.u[0].dof[i],
                                                                         self.u[1].dof[i],
                                                                         self.u[2].dof[i],
                                                                         self.u[3].dof[i],
                                                                         self.u[4].dof[i])
            if h is not None:
                self.u[0].dof[i] = h
            if hu is not None:
                self.u[1].dof[i] = hu
            if hv is not None:
                self.u[2].dof[i] = hv
            if heta is not None:
                self.u[3].dof[i] = heta
            if hw is not None:
                self.u[4].dof[i] = hw
    #

    def updateReflectingBoundaryConditions(self):
        self.forceStrongConditions = False
        for dummy, index in enumerate(self.boundaryIndex):
            vx = self.u[1].dof[index]
            vy = self.u[2].dof[index]
            vt = vx * self.normaly[index] - vy * self.normalx[index]
            self.u[1].dof[index] = vt * self.normaly[index]
            self.u[2].dof[index] = -vt * self.normalx[index]
    #

    def initDataStructures(self):
        # old vectors
        self.h_dof_old = numpy.copy(self.u[0].dof)
        self.hu_dof_old = numpy.copy(self.u[1].dof)
        self.hv_dof_old = numpy.copy(self.u[2].dof)
        self.heta_dof_old = numpy.copy(self.u[3].dof)
        self.hw_dof_old = numpy.copy(self.u[4].dof)
        # hEps
        eps = 1E-5
        self.hEps = eps * self.u[0].dof.max()
        # normal vectors
        self.normalx = numpy.zeros(self.u[0].dof.shape, 'd')
        self.normaly = numpy.zeros(self.u[0].dof.shape, 'd')
        # quantDOFs
        self.quantDOFs = numpy.zeros(self.u[0].dof.shape, 'd')
        # boundary Index: I do this in preStep since I need normalx and normaly to be initialized first
        # Allocate space for dLow (for the first stage in the SSP method)
        self.dLow = numpy.zeros(self.Cx.shape, 'd')
        self.hBT = numpy.zeros(self.Cx.shape, 'd')
        self.huBT = numpy.zeros(self.Cx.shape, 'd')
        self.hvBT = numpy.zeros(self.Cx.shape, 'd')
        self.hetaBT = numpy.zeros(self.Cx.shape, 'd')
        self.hwBT = numpy.zeros(self.Cx.shape, 'd')
        # get coordinates of DOFs
        self.getDOFsCoord()
        # some vectors for convex limiting
        self.dH_minus_dL = np.zeros(self.Cx.shape, 'd')
        self.muH_minus_muL = np.zeros(self.Cx.shape, 'd')
        self.extendedSourceTerm_hu = numpy.zeros(self.u[0].dof.shape, 'd')
        self.extendedSourceTerm_hv = numpy.zeros(self.u[0].dof.shape, 'd')
        self.extendedSourceTerm_heta = numpy.zeros(self.u[0].dof.shape, 'd')
        self.extendedSourceTerm_hw = numpy.zeros(self.u[0].dof.shape, 'd')
        self.new_SourceTerm_hu = numpy.zeros(self.u[0].dof.shape, 'd')
        self.new_SourceTerm_hv = numpy.zeros(self.u[0].dof.shape, 'd')
        self.new_SourceTerm_heta = numpy.zeros(self.u[0].dof.shape, 'd')
        self.new_SourceTerm_hw = numpy.zeros(self.u[0].dof.shape, 'd')
        self.dataStructuresInitialized = True
    #

    def getResidual(self, u, r):
        """
        Calculate the element residuals and add in to the global residual
        """

        # COMPUTE C MATRIX #
        if self.cterm_global is None:
            self.getCMatrices()

        # INIT DATA STRUCTURES #
        if self.dataStructuresInitialized == False:
            self.initDataStructures()

        # LOAD THE UNKNOWNS INTO THE FINITE ELEMENT DOFs #
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)

        # SET TO ZERO RESIDUAL #
        r.fill(0.0)

        # REFLECTING BOUNDARY CONDITIONS #
        if self.reflectingBoundaryConditions and self.boundaryIndex is not None:
            self.updateReflectingBoundaryConditions()
        #
        # INIT BOUNDARY INDEX #
        # NOTE: this must be done after the first call to getResidual (to have normalx and normaly initialized)
        # I do this in preStep if firstStep=True

        # CONSTRAINT DOFs #
        if self.coefficients.constrainedDOFs is not None:
            self.updateConstrainedDOFs()
        #
        # DIRICHLET BOUNDARY CONDITIONS #
        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.items()):
                    self.u[cj].dof[dofN] = g(
                        self.dirichletConditionsForceDOF[cj].DOFBoundaryPointDict[dofN], self.timeIntegration.t)
        #
        # CHECK POSITIVITY OF WATER HEIGHT # changed to 1E-4 -EJT
        if (self.check_positivity_water_height == False):
            assert self.u[0].dof.min(
            ) > -1E-5 * self.u[0].dof.max(), ("Negative water height: ", self.u[0].dof.min())

        self.calculateResidual(
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
            self.u[0].femSpace.dofMap.l2g,
            self.u[1].femSpace.dofMap.l2g,
            self.h_dof_old,
            self.hu_dof_old,
            self.hv_dof_old,
            self.heta_dof_old,
            self.hw_dof_old,
            self.coefficients.b.dof,
            self.u[0].dof,
            self.u[1].dof,
            self.u[2].dof,
            self.u[3].dof,
            self.u[4].dof,
            self.h_dof_sge,
            self.hu_dof_sge,
            self.hv_dof_sge,
            self.heta_dof_sge,
            self.hw_dof_sge,
            self.timeIntegration.m_tmp[0],
            self.timeIntegration.m_tmp[1],
            self.timeIntegration.m_tmp[2],
            self.q[('f', 0)],
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.q[('cfl', 0)],
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
            self.offset[3],
            self.offset[4],
            self.stride[0],
            self.stride[1],
            self.stride[2],
            self.stride[3],
            self.stride[4],
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
            self.elementResidual[0],
            self.Cx,
            self.Cy,
            self.CTx,
            self.CTy,
            self.numDOFsPerEqn,
            self.numNonZeroEntries,
            self.rowptr_cMatrix,
            self.colind_cMatrix,
            self.ML,
            self.timeIntegration.runCFL,
            self.hEps,
            self.hReg,
            self.q[('u', 0)],
            self.q[('u', 1)],
            self.q[('u', 2)],
            self.q[('u', 3)],
            self.q[('u', 4)],
            self.extendedSourceTerm_hu,
            self.extendedSourceTerm_hv,
            self.extendedSourceTerm_heta,
            self.extendedSourceTerm_hw,
            self.dH_minus_dL,
            self.muH_minus_muL,
            self.coefficients.cE,
            self.coefficients.LUMPED_MASS_MATRIX,
            self.timeIntegration.dt,
            self.coefficients.LINEAR_FRICTION,
            self.coefficients.mannings,
            self.quantDOFs,
            self.secondCallCalculateResidual,
            self.COMPUTE_NORMALS,
            self.normalx,
            self.normaly,
            self.dLow,
            self.hBT,
            self.huBT,
            self.hvBT,
            self.hetaBT,
            self.hwBT,
            self.timeIntegration.lstage,
            self.new_SourceTerm_hu,
            self.new_SourceTerm_hv,
            self.new_SourceTerm_heta,
            self.new_SourceTerm_hw)

        self.COMPUTE_NORMALS = 0
        if self.forceStrongConditions:
            for cj in range(len(self.dirichletConditionsForceDOF)):
                for dofN, g in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.items()):
                    r[self.offset[cj] + self.stride[cj] * dofN] = 0.
        #
        if self.constrainedDOFsIndices is not None:
            for index in self.constrainedDOFsIndices:
                for cj in range(self.nc):
                    global_dofN = self.offset[cj] + self.stride[cj] * index
                    r[global_dofN] = 0.
        #
        logEvent("Global residual", level=9, data=r)
        # mwf decide if this is reasonable for keeping solver statistics
        self.nonlinear_function_evaluations += 1

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        self.calculateJacobian(
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
            self.u[3].dof,
            self.u[4].dof,
            self.h_dof_sge,
            self.hu_dof_sge,
            self.hv_dof_sge,
            self.heta_dof_sge,
            self.hw_dof_sge,
            self.timeIntegration.beta_bdf[0],
            self.timeIntegration.beta_bdf[1],
            self.timeIntegration.beta_bdf[2],
            self.q[('cfl', 0)],
            self.coefficients.sdInfo[(1, 1)][0],
            self.coefficients.sdInfo[(1, 1)][1],
            self.coefficients.sdInfo[(1, 2)][0],
            self.coefficients.sdInfo[(1, 2)][1],
            self.coefficients.sdInfo[(2, 2)][0],
            self.coefficients.sdInfo[(2, 2)][1],
            self.coefficients.sdInfo[(2, 1)][0],
            self.coefficients.sdInfo[(2, 1)][1],
            # h
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            self.csrRowIndeces[(0, 1)], self.csrColumnOffsets[(0, 1)],
            self.csrRowIndeces[(0, 2)], self.csrColumnOffsets[(0, 2)],
            self.csrRowIndeces[(0, 3)], self.csrColumnOffsets[(0, 3)],
            self.csrRowIndeces[(0, 4)], self.csrColumnOffsets[(0, 4)],
            # hu
            self.csrRowIndeces[(1, 0)], self.csrColumnOffsets[(1, 0)],
            self.csrRowIndeces[(1, 1)], self.csrColumnOffsets[(1, 1)],
            self.csrRowIndeces[(1, 2)], self.csrColumnOffsets[(1, 2)],
            self.csrRowIndeces[(1, 3)], self.csrColumnOffsets[(1, 3)],
            self.csrRowIndeces[(1, 4)], self.csrColumnOffsets[(1, 4)],
            # hv
            self.csrRowIndeces[(2, 0)], self.csrColumnOffsets[(2, 0)],
            self.csrRowIndeces[(2, 1)], self.csrColumnOffsets[(2, 1)],
            self.csrRowIndeces[(2, 2)], self.csrColumnOffsets[(2, 2)],
            self.csrRowIndeces[(2, 3)], self.csrColumnOffsets[(2, 3)],
            self.csrRowIndeces[(2, 4)], self.csrColumnOffsets[(2, 4)],
            # heta
            self.csrRowIndeces[(3, 0)], self.csrColumnOffsets[(3, 0)],
            self.csrRowIndeces[(3, 1)], self.csrColumnOffsets[(3, 1)],
            self.csrRowIndeces[(3, 2)], self.csrColumnOffsets[(3, 2)],
            self.csrRowIndeces[(3, 3)], self.csrColumnOffsets[(3, 3)],
            self.csrRowIndeces[(3, 4)], self.csrColumnOffsets[(3, 4)],
            # hw
            self.csrRowIndeces[(4, 0)], self.csrColumnOffsets[(4, 0)],
            self.csrRowIndeces[(4, 1)], self.csrColumnOffsets[(4, 1)],
            self.csrRowIndeces[(4, 2)], self.csrColumnOffsets[(4, 2)],
            self.csrRowIndeces[(4, 3)], self.csrColumnOffsets[(4, 3)],
            self.csrRowIndeces[(4, 4)], self.csrColumnOffsets[(4, 4)],
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
            self.csrColumnOffsets_eb[(2, 2)],
            self.timeIntegration.dt)

        # Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            scaling = 1.0  # probably want to add some scaling to match non-dirichlet diagonals in linear system
            for cj in range(self.nc):
                for dofN in list(self.dirichletConditionsForceDOF[cj].DOFBoundaryConditionsDict.keys()):
                    global_dofN = self.offset[cj] + self.stride[cj] * dofN
                    for i in range(self.rowptr[global_dofN], self.rowptr[global_dofN + 1]):
                        if (self.colind[i] == global_dofN):
                            self.nzval[i] = scaling
                        else:
                            self.nzval[i] = 0.0
        #
        if self.constrainedDOFsIndices is not None:
            for index in self.constrainedDOFsIndices:
                for cj in range(self.nc):
                    global_dofN = self.offset[cj] + self.stride[cj] * index
                    for i in range(self.rowptr[global_dofN], self.rowptr[global_dofN + 1]):
                        if (self.colind[i] == global_dofN):
                            self.nzval[i] = 1.0
                        else:
                            self.nzval[i] = 0.0
        #
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
        self.u[0].femSpace.elementMaps.getValues(
            self.elementQuadraturePoints, self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(
            self.timeIntegration.t, self.q)

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
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[1].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        self.fluxBoundaryConditionsObjectsDict = dict([(cj, FluxBoundaryConditions(self.mesh,
                                                                                   self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                   self.ebqe[(
                                                                                       'x')],
                                                                                   self.advectiveFluxBoundaryConditionsSetterDict[
                                                                                       cj],
                                                                                   self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in list(self.advectiveFluxBoundaryConditionsSetterDict.keys())])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(
            self.timeIntegration.t, self.ebqe)

    def estimate_mt(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        self.h_dof_sge[:] = self.u[0].dof
        self.hu_dof_sge[:] = self.u[1].dof
        self.hv_dof_sge[:] = self.u[2].dof
        self.heta_dof_sge[:] = self.u[3].dof
        self.hw_dof_sge[:] = self.u[4].dof
        OneLevelTransport.calculateAuxiliaryQuantitiesAfterStep(self)

    def getForce(self, cg, forceExtractionFaces, force, moment):
        pass
