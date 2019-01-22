from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import proteus
from proteus.mprans.cBlendedSpaces import *


class NumericalFlux(proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior):
    def __init__(self, vt, getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions,
                 getPeriodicBoundaryConditions=None):
        proteus.NumericalFlux.Advection_DiagonalUpwind_Diffusion_IIPG_exterior.__init__(self, vt,
                                                                                        getPointwiseBoundaryConditions,
                                                                                        getAdvectiveFluxBoundaryConditions,
                                                                                        getDiffusiveFluxBoundaryConditions)

class Coefficients(proteus.TransportCoefficients.TC_base):
    from proteus.ctransportCoefficients import VOFCoefficientsEvaluate
    from proteus.ctransportCoefficients import VolumeAveragedVOFCoefficientsEvaluate
    from proteus.cfemIntegrals import copyExteriorElementBoundaryValuesFromElementBoundaryValues

    def __init__(self,
                 he=0.1,
                 ME_model=0,
                 epsFact=0.0,
                 forceStrongConditions=0,
                 # OUTPUT quantDOFs
                 outputQuantDOFs=False,
                 #NULLSPACE INFO
                 epsilon=0.01,
                 ALPHA_FOR_GALERKIN_SOLUTION=1,
                 AUTOMATED_ALPHA=True,
                 PROBLEM_TYPE=0):

        self.ALPHA_FOR_GALERKIN_SOLUTION=ALPHA_FOR_GALERKIN_SOLUTION
        self.AUTOMATED_ALPHA=AUTOMATED_ALPHA
        self.he=he
        self.PROBLEM_TYPE=PROBLEM_TYPE

        self.DISCRETE_UPWINDING = 1
        if self.PROBLEM_TYPE==0 and self.AUTOMATED_ALPHA==0:
            self.DISCRETE_UPWINDING = 0

        self.epsilon=epsilon
        self.variableNames = ['u']
        nc = 1
        mass = {0: {0: 'linear'}}
        advection = {0: {0: 'linear'}}
        hamiltonian = {}
        diffusion = {}
        potential = {}
        reaction = {}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian,
                         self.variableNames)
        self.epsFact = epsFact
        self.modelIndex = ME_model
        self.forceStrongConditions = forceStrongConditions
        self.outputQuantDOFs = outputQuantDOFs

    def initializeMesh(self, mesh):
        self.eps = self.epsFact * mesh.h

    def attachModels(self, modelList):
        self.model = modelList[self.modelIndex]
        self.q_v = np.zeros(self.model.q[('grad(u)',0)].shape,'d')
        self.ebqe_v = np.zeros(self.model.ebqe[('grad(u)',0)].shape,'d')
        self.ebqe_phi = np.zeros(self.model.ebqe[('u',0)].shape,'d')#cek hack, we don't need this

    def initializeElementQuadrature(self, t, cq):
        pass

    def initializeElementBoundaryQuadrature(self, t, cebq, cebq_global):
        pass

    def initializeGlobalExteriorElementBoundaryQuadrature(self, t, cebqe):
        pass

    def preStep(self, t, firstStep=False):
        # SAVE OLD SOLUTION #
        self.model.u_dof_old[:] = self.model.u[0].dof

        if self.model.hasForceFieldAsFunction:
            self.model.updateForceFieldAsFunction()

        # COMPUTE BLENDING FUNCTION (if given by user)
        if self.model.hasBlendingFunction:
            self.model.updateBlendingFunction()

        # COMPUTE NEW VELOCITY (if given by user) #
        if self.model.hasVelocityFieldAsFunction:
            self.model.updateVelocityFieldAsFunction()

        # COMPUTE FORCE FIELD (if given by user)
        if self.model.hasForceFieldAsFunction:
            self.model.updateForceFieldAsFunction()

        copyInstructions = {}
        return copyInstructions

    def postStep(self, t, firstStep=False):
        print "********************... Min, Max: ", self.model.u[0].dof.min(), self.model.u[0].dof.max()
        print "********************... Ndof: ", len(self.model.u[0].dof)
        tolerance = 1.0E-10
        MIN=-1.0
        MAX=1.0
        #self.model.quantDOFs[:] = 1.0*(self.model.u[0].dof > MAX+tolerance) - 1.0*(self.model.u[0].dof < MIN-tolerance)
        #self.model.quantDOFs[:] = 1.0*(self.model.u[0].dof < MIN-tolerance)
        #self.model.quantDOFs[:] = 1.0*(self.model.u[0].dof >= MAX+tolerance)

        # COMPUTE NORMS #
        self.model.getMetricsAtEOS()
        self.model.metricsAtEOS.write(repr(np.sqrt(self.model.global_L2))+","+
                                      repr(np.sqrt(self.model.global_H1))+","+
                                      repr(np.sqrt(self.model.global_L2_Omega1))+","+
                                      repr(np.sqrt(self.model.global_H1_Omega1))+","+
                                      repr(np.sqrt(self.model.global_L2_Omega2))+","+
                                      repr(np.sqrt(self.model.global_H1_Omega2))+","+
                                      repr(np.sqrt(self.model.global_L2_sH))+","+
                                      repr(np.sqrt(self.model.global_L2_1msH))+"\n")
        self.model.metricsAtEOS.flush()

        copyInstructions = {}
        return copyInstructions

    def updateToMovingDomain(self, t, c):
        pass

    def evaluate(self, t, c):
        pass

class LevelModel(proteus.Transport.OneLevelTransport):
    nCalls = 0

    def __init__(self,
                 uDict,
                 phiDict,
                 auxDict,
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
                 stressTraceBoundaryConditionsSetterDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='defaultName',
                 reuse_trial_and_test_quadrature=True,
                 sd=True,
                 movingDomain=False,
                 bdyNullSpace=False):

        self.auxiliaryCallCalculateResidual = False
        #
        # set the objects describing the method and boundary conditions
        #
        self.bdyNullSpace = bdyNullSpace
        self.movingDomain = movingDomain
        self.tLast_mesh = None
        #
        self.name = name
        self.sd = sd
        self.Hess = False
        self.lowmem = True
        self.timeTerm = True  # allow turning off  the  time derivative
        # self.lowmem=False
        self.testIsTrial = True
        self.phiTrialIsTrial = True
        self.u = uDict
        self.aux = auxDict
        self.ua = {}  # analytical solutions
        self.phi = phiDict
        self.dphi = {}
        self.matType = matType
        # mwf try to reuse test and trial information across components if spaces are the same
        self.reuse_test_trial_quadrature = reuse_trial_and_test_quadrature  # True#False
        if self.reuse_test_trial_quadrature:
            for ci in range(1, coefficients.nc):
                assert self.u[ci].femSpace.__class__.__name__ == self.u[0].femSpace.__class__.__name__, "to reuse_test_trial_quad all femSpaces must be the same!"
        self.u_dof_old = None

        # Simplicial Mesh
        self.mesh = self.u[0].femSpace.mesh  # assume the same mesh for  all components for now
        self.testSpace = testSpaceDict
        self.dirichletConditions = dofBoundaryConditionsDict
        self.dirichletNodeSetList = None  # explicit Dirichlet  conditions for now, no Dirichlet BC constraints
        self.coefficients = coefficients
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing


        self.conservativeFlux = conservativeFluxDict  # no velocity post-processing for now
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        if self.stabilization is not None:
            raise NotImplementedError ("No stabilization for this model")

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
            raise NotImplementedError ("No stabilization for this model")

        if self.shockCapturing is not None:
            raise NotImplementedError ("No shock capturing for this model")

        if massLumping:
            raise NotImplementedError ("N/A for this model")

        if reactionLumping:
            raise NotImplementedError ("N/A for this model")

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
        #
        # storage dictionaries
        #
        self.scalars_element = set()
        #
        # simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q = {}
        self.ebq = {}
        self.ebq_global = {}
        self.ebqe = {}
        self.phi_ip = {}
        # mesh
        self.q['x'] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, 3), 'd')
        self.ebqe['x'] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary, 3), 'd')
        self.q[('u', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('dV_u', 0)] = (old_div(1.0, self.mesh.nElements_global)) * numpy.ones((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('grad(u)', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element, self.nSpace_global), 'd')
        self.q[('m_last', 0)] = numpy.zeros((self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[('m_tmp', 0)] = self.q[('u', 0)]
        self.ebqe[('u', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('grad(u)', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                 self.nElementBoundaryQuadraturePoints_elementBoundary, self.nSpace_global), 'd')
        self.ebqe[('advectiveFlux_bc_flag', 0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'i')
        self.ebqe[('advectiveFlux_bc', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
        self.ebqe[('advectiveFlux', 0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global, self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')

        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set([('u', ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature = set()
        self.tensors_elementBoundaryQuadrature = set()
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
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(self)
        logEvent(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        # TODO get rid of this
        for ci, fbcObject in list(self.fluxBoundaryConditionsObjectsDict.items()):
            self.ebqe[('advectiveFlux_bc_flag', ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc', ci)].shape, 'i')
            for t, g in list(fbcObject.advectiveFluxBoundaryConditionsDict.items()):
                if ci in self.coefficients.advection:
                    self.ebqe[('advectiveFlux_bc', ci)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag', ci)][t[0], t[1]] = 1

        if hasattr(self.numericalFlux, 'setDirichletValues'):
            self.numericalFlux.setDirichletValues(self.ebqe)
        if not hasattr(self.numericalFlux, 'isDOFBoundary'):
            self.numericalFlux.isDOFBoundary = {0: numpy.zeros(self.ebqe[('u', 0)].shape, 'i')}
        if not hasattr(self.numericalFlux, 'ebqe'):
            self.numericalFlux.ebqe = {('u', 0): numpy.zeros(self.ebqe[('u', 0)].shape, 'd')}
        # TODO how to handle redistancing calls for calculateCoefficients,calculateElementResidual etc
        self.globalResidualDummy = None
        compKernelFlag = 0
        self.blendedSpaces = cBlendedSpaces_base(self.nSpace_global,
                                                 self.nQuadraturePoints_element,
                                                 self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
                                                 self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
                                                 self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
                                                 self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                 compKernelFlag)

        self.forceStrongConditions=self.coefficients.forceStrongConditions
        if self.forceStrongConditions:
            self.dirichletConditionsForceDOF = DOFBoundaryConditions(self.u[0].femSpace, dofBoundaryConditionsSetterDict[0], weakDirichletConditions=False)
        #

        # mql. Allow the user to provide functions to define the velocity field
        #Velocity field
        self.hasVelocityFieldAsFunction = False
        if ('velocityFieldAsFunction') in dir(options):
            self.velocityFieldAsFunction = options.velocityFieldAsFunction
            self.hasVelocityFieldAsFunction = True

        #Blending function
        self.hasBlendingFunction = False
        if ('blendingFunction') in dir (options):
            self.blendingFunction = options.blendingFunction
            self.hasBlendingFunction = True

        #Exact solution
        self.hasAlphaFunction = False
        if ('alphaFunction') in dir(options):
            self.alphaFunction = options.alphaFunction
            self.hasAlphaFunction = True
        self.q['alpha_value'] = numpy.zeros((self.mesh.nElements_global,
                                             self.nQuadraturePoints_element),'d')
        #Force field
        self.hasForceFieldAsFunction = False
        if ('forceFieldAsFunction') in dir(options):
            self.forceFieldAsFunction = options.forceFieldAsFunction
            self.hasForceFieldAsFunction = True
        self.q['force'] = numpy.zeros((self.mesh.nElements_global,
                                           self.nQuadraturePoints_element),'d')

        # mql. For edge based stabilization
        self.is_dof_external = numpy.zeros(self.u[0].dof.shape, 'd')
        self.is_dof_internal = numpy.zeros(self.u[0].dof.shape, 'd')

        self.is_boundary = numpy.zeros(self.u[0].dof.shape, 'd')
        self.is_cell_boundary = numpy.zeros(self.u[0].dof.shape, 'd')
        self.num_hi = numpy.zeros(self.u[0].dof.shape, 'd')
        self.den_hi = numpy.zeros(self.u[0].dof.shape, 'd')
        self.GALERKIN_SOLUTION=0
        self.dLow=None
        self.quantDOFs = numpy.zeros(self.u[0].dof.shape, 'd')
        self.quantDOFs2 = numpy.zeros(self.u[0].dof.shape, 'd')
        self.quantDOFs3 = numpy.zeros(self.u[0].dof.shape, 'd')
        self.galerkin_solution = numpy.zeros(self.u[0].dof.shape, 'd')
        self.boundaryValues = None
        self.isBoundary = None
        self.hi = numpy.zeros(self.u[0].dof.shape, 'd')
        self.gamma_dof = numpy.zeros(self.u[0].dof.shape, 'd')
        self.beta_dof = numpy.zeros(self.u[0].dof.shape, 'd')
        self.element_He = numpy.zeros(self.mesh.nElements_global, 'd')

        # mql. For blending functions
        self.blendingFunctionDOFs = numpy.zeros(self.u[0].dof.shape,'d')
        # get Q2mesh_dof. This has the location of all DOFs in the Q2 mesh
        self.xCoord_dof = numpy.zeros(self.u[0].dof.shape,'d')
        self.yCoord_dof = numpy.zeros(self.u[0].dof.shape,'d')
        self.zCoord_dof = numpy.zeros(self.u[0].dof.shape,'d')
        for eN in range(self.mesh.nElements_global):
            for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                gi = self.offset[0]+self.stride[0]*self.u[0].femSpace.dofMap.l2g[eN,i]
                self.xCoord_dof[gi] = self.u[0].femSpace.interpolationPoints[eN,i,0]
                self.yCoord_dof[gi] = self.u[0].femSpace.interpolationPoints[eN,i,1]
                self.zCoord_dof[gi] = self.u[0].femSpace.interpolationPoints[eN,i,2]
        #
        # mql. METRICS #
        self.exactSolution = options.exactSolution
        self.exactGrad = options.exactGrad
        self.global_L2 = 0.0
        self.global_H1 = 0.0
        self.global_L2_Omega1 = 0.0
        self.global_H1_Omega1 = 0.0
        self.global_L2_Omega2 = 0.0
        self.global_H1_Omega2 = 0.0
        self.global_L2_sH = 0.0
        self.global_L2_1msH = 0.0
        self.metricsAtEOS = open(self.name+"_metricsAtEOS.csv","w")
        self.metricsAtEOS.write('global_L2'+","+
                                'global_H1'+","+
                                'global_L2_Omega1'+","+
                                'global_H1_Omega1'+","+
                                'global_L2_Omega2'+","+
                                'global_H1_Omega2'+","+
                                'global_L2_sH'+","
                                'global_L2_1msH'+"\n")

    def getSmoothnessIndicator(self,u):
        self.num_hi[:]=0
        self.den_hi[:]=0
        self.is_dof_external[:]=0
        self.is_dof_internal[:]=0

        self.blendedSpaces.calculateSmoothnessIndicator(  # element
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.Hessian_psi,
            # physics
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.l2g[0]['freeGlobal'],
            u,
            self.offset[0], self.stride[0],
            self.is_dof_external,
            self.is_dof_internal,
            self.num_hi,
            self.den_hi,
            self.hi,
            self.element_He,
            self.coefficients.he,
            self.xCoord_dof,
            self.yCoord_dof,
            self.zCoord_dof,
            self.gamma_dof,
            self.beta_dof,
            len(self.rowptr)-1, # num of DOFs
            self.nnz,
            self.rowptr,
            self.colind)

    def getMetricsAtEOS(self):
        """
        Calculate some metrics at EOS (End Of Simulation)
        """
        u_exact = numpy.zeros(self.q[('u',0)].shape,'d')
        gradx_u_exact = numpy.zeros(self.q[('u',0)].shape,'d')
        grady_u_exact = numpy.zeros(self.q[('u',0)].shape,'d')
        X = {0:self.q[('x')][:,:,0],
             1:self.q[('x')][:,:,1],
             2:self.q[('x')][:,:,2]}
        u_exact[:] = self.exactSolution[0](X,0)
        gradx_u_exact[:] = self.exactGrad[0](X,0)
        grady_u_exact[:] = self.exactGrad[1](X,0)

        meshSize = np.zeros(1)
        meshSize[0] = 1.0/(np.sqrt(self.u[0].dof.size)-1)
        (global_L2,
         global_H1,
         global_L2_Omega1,
         global_H1_Omega1,
         global_Omega1,
         global_L2_Omega2,
         global_H1_Omega2,
         global_Omega2,
         global_L2_sH,
         global_L2_1msH) = self.blendedSpaces.calculateMetricsAtEOS(#element
             self.u[0].femSpace.elementMaps.psi,
             self.u[0].femSpace.elementMaps.grad_psi,
             self.mesh.nodeArray,
             self.mesh.elementNodesArray,
             self.elementQuadratureWeights[('u',0)],
             self.u[0].femSpace.psi,
             self.u[0].femSpace.grad_psi,
             self.u[0].femSpace.psi,
             #physics
             self.mesh.nElements_global,
             self.mesh.nElements_owned,
             self.u[0].femSpace.dofMap.l2g,
             self.mesh.elementDiametersArray,
             meshSize,
             self.mesh.nodeDiametersArray,
             2.0, #epsFactHeaviside
             self.q[('u',0)],
             self.q[('grad(u)',0)],
             u_exact,
             gradx_u_exact,
             grady_u_exact,
             self.offset[0],self.stride[0])


        from proteus.flcbdfWrappers import globalSum
        # Interface metrics
        relative_norms=False
        self.global_L2 = globalSum(global_L2)
        self.global_H1 = globalSum(global_H1)
        self.global_L2_sH = globalSum(global_L2_sH)
        if relative_norms:
        #
            self.global_L2_Omega1 = globalSum(global_L2_Omega1)/globalSum(global_Omega1)
            self.global_H1_Omega1 = globalSum(global_H1_Omega1)/globalSum(global_Omega1)
            self.global_L2_Omega2 = globalSum(global_L2_Omega2)/globalSum(global_Omega2)
            self.global_H1_Omega2 = globalSum(global_H1_Omega2)/globalSum(global_Omega2)
            self.global_L2_1msH = globalSum(global_L2_sH)/globalSum(global_L2_1msH)
        #
        else:
            self.global_L2_Omega1 = globalSum(global_L2_Omega1)
            self.global_H1_Omega1 = globalSum(global_H1_Omega1)
            self.global_L2_Omega2 = globalSum(global_L2_Omega2)
            self.global_H1_Omega2 = globalSum(global_H1_Omega2)
            self.global_L2_1msH = globalSum(global_L2_1msH)
        #

    def calculateCoefficients(self):
        pass

    def updateBlendingFunction(self):
        # DO LUMPED L2 PROJECTION #
        self.blendingFunctionViaLumpedL2Projection=False
        if self.blendingFunctionViaLumpedL2Projection:
            logEvent("Computing blending function via lumped L2 projection", level=2)
            self.updateAlphaFunction()
            self.blendedSpaces.getLumpedL2Projection(#element
                self.u[0].femSpace.elementMaps.psi,
                self.u[0].femSpace.elementMaps.grad_psi,
                self.mesh.nodeArray,
                self.mesh.elementNodesArray,
                self.elementQuadratureWeights[('u',0)],
                self.aux[0].femSpace.psi,
                self.aux[0].femSpace.grad_psi,
                self.aux[0].femSpace.psi,
                self.mesh.nElements_global,
                self.u[0].femSpace.dofMap.l2g,
                self.mesh.elementDiametersArray,
                self.q['alpha_value'],
                self.offset[0],self.stride[0],
                self.nFreeDOF_global[0], #numDOFs
                self.blendingFunctionDOFs)
        else:
            # DO FINITE ELEMENT INTERPOLATION #
            logEvent("Computing blending function via FE interpolation", level=2)
            self.x = numpy.zeros(self.u[0].dof.shape,'d')
            self.y = numpy.zeros(self.u[0].dof.shape,'d')
            self.z = numpy.zeros(self.u[0].dof.shape,'d')
            for eN in range(self.mesh.nElements_global):
                for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                    gi = self.offset[0]+self.stride[0]*self.u[0].femSpace.dofMap.l2g[eN,i]
                    self.x[gi] = self.u[0].femSpace.interpolationPoints[eN,i,0]
                    self.y[gi] = self.u[0].femSpace.interpolationPoints[eN,i,1]
                    self.z[gi] = self.u[0].femSpace.interpolationPoints[eN,i,2]

            X = {0:self.x,
                 1:self.y,
                 2:self.z}
            t = self.timeIntegration.t

            self.blendingFunctionDOFs[:] = self.blendingFunction[0](X,t)

    def updateVelocityFieldAsFunction(self):
        X = {0: self.q[('x')][:, :, 0],
             1: self.q[('x')][:, :, 1],
             2: self.q[('x')][:, :, 2]}
        t = self.timeIntegration.t
        self.coefficients.q_v[..., 0] = self.velocityFieldAsFunction[0](X, t)
        self.coefficients.q_v[..., 1] = self.velocityFieldAsFunction[1](X, t)
        if (self.nSpace_global == 3):
            self.coefficients.q_v[..., 2] = self.velocityFieldAsFunction[2](X, t)

        # BOUNDARY
        ebqe_X = {0: self.ebqe['x'][:, :, 0],
                  1: self.ebqe['x'][:, :, 1],
                  2: self.ebqe['x'][:, :, 2]}
        self.coefficients.ebqe_v[..., 0] = self.velocityFieldAsFunction[0](ebqe_X, t)
        self.coefficients.ebqe_v[..., 1] = self.velocityFieldAsFunction[1](ebqe_X, t)
        if (self.nSpace_global == 3):
            self.coefficients.ebqe_v[..., 2] = self.velocityFieldAsFunction[2](ebqe_X, t)

    def updateForceFieldAsFunction(self):
        X = {0:self.q[('x')][:,:,0],
             1:self.q[('x')][:,:,1],
             2:self.q[('x')][:,:,2]}
        t = self.timeIntegration.t
        self.q['force'][:] = self.forceFieldAsFunction[0](X,t)

    def updateAlphaFunction(self):
        X = {0:self.q[('x')][:,:,0],
             1:self.q[('x')][:,:,1],
             2:self.q[('x')][:,:,2]}
        t = self.timeIntegration.t
        self.q['alpha_value'][:] = self.alphaFunction[0](X,t)

    def calculateElementResidual(self):
        if self.globalResidualDummy is not None:
            self.getResidual(self.u[0].dof, self.globalResidualDummy)

    def getBoundaryValues(self):
        # get x,y coordinates of all DOFs #
        l2g = self.l2g[0]['freeGlobal']
        self.dofsXCoord = numpy.zeros(self.u[0].dof.shape,'d')
        self.dofsYCoord = numpy.zeros(self.u[0].dof.shape,'d')
        for eN in range(self.mesh.nElements_global):
            for k in range(self.nQuadraturePoints_element):
                for i in range(self.nDOF_test_element[0]):
                    gi = self.offset[0]+self.stride[0]*l2g[eN,i]
                    self.dofsXCoord[gi] = self.u[0].femSpace.interpolationPoints[eN,i,0]
                    self.dofsYCoord[gi] = self.u[0].femSpace.interpolationPoints[eN,i,1]
        #
        # Get DOFs values #
        self.isBoundary = numpy.zeros(self.u[0].dof.shape,'d')
        self.boundaryValues = {} #numpy.zeros(0,'d')

        eps = 0.1/(np.sqrt(self.u[0].dof.size)-1)
        for gi in range(self.u[0].dof.shape[0]):
            x = self.dofsXCoord[gi]
            y = self.dofsYCoord[gi]
            # OUTTER BOUNDARY #
            if x==0.0 or x==1.0 or y==0.0 or y==1.0:
                value=-1.0
                self.boundaryValues[gi] = value
                self.isBoundary[gi] = 1.0
            # INTERNAL BOUNDARY #
            if (x >= 4./9-eps and x <= 5./9+eps and y >= 4./9-eps and y <= 5./9+eps):
                value=1.0
                self.boundaryValues[gi] = value
                self.isBoundary[gi] = 1.0
            #if x==0:
            #    self.boundaryValues[gi] = -1.
            #    self.isBoundary[gi]=1.0
            #elif x==1:
            #    self.boundaryValues[gi] = 1.
            #    self.isBoundary[gi]=1.0
        #

    def getResidual(self, u, r):
        import pdb
        import copy
        """
        Calculate the element residuals and add in to the global residual
        """
        if self.boundaryValues is None:
            self.getBoundaryValues()
        # JACOBIANS (FOR ELEMENT TRANSFORMATION)
        self.q[('J')] = np.zeros((self.mesh.nElements_global,
                                  self.nQuadraturePoints_element,
                                  self.nSpace_global,
                                  self.nSpace_global),
                                 'd')
        self.q[('inverse(J)')] = np.zeros((self.mesh.nElements_global,
                                           self.nQuadraturePoints_element,
                                           self.nSpace_global,
                                           self.nSpace_global),
                                          'd')
        self.q[('det(J)')] = np.zeros((self.mesh.nElements_global,
                                       self.nQuadraturePoints_element),
                                      'd')
        self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                         self.q['J'],
                                                         self.q['inverse(J)'],
                                                         self.q['det(J)'])
        self.q['abs(det(J))'] = np.abs(self.q['det(J)'])

        if self.u_dof_old is None:
            # Pass initial condition to u_dof_old
            self.u_dof_old = numpy.copy(self.u[0].dof)

        # This is dummy. I just care about the csr structure of the sparse matrix
        self.dLow = np.zeros(self.nnz,'d')

        r.fill(0.0)
        # Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)
        # cek can put in logic to skip of BC's don't depend on t or u
        # Dirichlet boundary conditions
        # if hasattr(self.numericalFlux,'setDirichletValues'):
        self.numericalFlux.setDirichletValues(self.ebqe)
        # flux boundary conditions
        for t, g in list(self.fluxBoundaryConditionsObjectsDict[0].advectiveFluxBoundaryConditionsDict.items()):
            self.ebqe[('advectiveFlux_bc', 0)][t[0], t[1]] = g(self.ebqe[('x')][t[0], t[1]], self.timeIntegration.t)
            self.ebqe[('advectiveFlux_bc_flag', 0)][t[0], t[1]] = 1

        if self.forceStrongConditions:
            for dofN, g in list(self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.items()):
                self.u[0].dof[dofN] = g(self.dirichletConditionsForceDOF.DOFBoundaryPointDict[dofN], self.timeIntegration.t)
                self.is_boundary[dofN] = 1
        #
        for eN in range(self.mesh.nElements_global):
            cell_has_boundary =False
            for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                gi = self.offset[0]+self.stride[0]*self.u[0].femSpace.dofMap.l2g[eN,i]
                if self.is_boundary[gi] == 1:
                    cell_has_boundary = True
            if cell_has_boundary:
                for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                    gi = self.offset[0]+self.stride[0]*self.u[0].femSpace.dofMap.l2g[eN,i]
                    self.is_cell_boundary[gi] = 1
        #
        self.calculateResidual = self.blendedSpaces.calculateResidual
        self.calculateJacobian = self.blendedSpaces.calculateJacobian

        self.dLow.fill(0.0)
        #import pdb; pdb.set_trace()
        #print (self.u[0].femSpace.Hessian_psi.min(), self.u[0].femSpace.Hessian_psi.max())
        self.calculateResidual(  # element
            self.timeIntegration.dt,
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            # physics
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.l2g[0]['freeGlobal'],
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.u_dof_old,  # For Backward Euler this is un, for SSP this is the lstage
            self.coefficients.q_v,
            self.q[('u', 0)],
            self.q[('grad(u)',0)],
            self.offset[0], self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.coefficients.ebqe_phi,
            self.coefficients.epsFact,
            self.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux', 0)],
            # ENTROPY VISCOSITY and ARTIFICIAL COMRPESSION
            # PARAMETERS FOR EDGE VISCOSITY
            len(self.rowptr)-1, # num of DOFs
            self.nnz,
            self.rowptr,  # Row indices for Sparsity Pattern (convenient for DOF loops)
            self.colind,  # Column indices for Sparsity Pattern (convenient for DOF loops)
            self.csrRowIndeces[(0, 0)],  # row indices (convenient for element loops)
            self.csrColumnOffsets[(0, 0)],  # column indices (convenient for element loops)
            self.csrColumnOffsets_eb[(0, 0)],  # indices for boundary terms
            # For blending spaces
            self.q['force'],
            self.blendingFunctionDOFs, #alpha
            self.aux[0].femSpace.psi,
            self.aux[0].femSpace.grad_psi,
            self.dLow,
            self.coefficients.epsilon,
            self.coefficients.AUTOMATED_ALPHA,
            self.coefficients.PROBLEM_TYPE,
            self.coefficients.DISCRETE_UPWINDING,
            self.coefficients.ALPHA_FOR_GALERKIN_SOLUTION,
            self.GALERKIN_SOLUTION,
            self.galerkin_solution,
            self.gamma_dof,
            self.beta_dof,
            self.quantDOFs,
            self.quantDOFs2,
            self.quantDOFs3,
            self.is_boundary,
            self.is_cell_boundary)

        if self.forceStrongConditions:
            for dofN, g in list(self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.items()):
                r[dofN] = 0

        logEvent("Global residual", level=9, data=r)
        self.nonlinear_function_evaluations += 1
        if self.globalResidualDummy is None:
            self.globalResidualDummy = numpy.zeros(r.shape, 'd')

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian,
                                       jacobian)
        self.calculateJacobian(  # element
            self.timeIntegration.dt,
            self.u[0].femSpace.elementMaps.psi,
            self.u[0].femSpace.elementMaps.grad_psi,
            self.mesh.nodeArray,
            self.mesh.elementNodesArray,
            self.elementQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            self.u[0].femSpace.psi,
            self.u[0].femSpace.grad_psi,
            # element boundary
            self.u[0].femSpace.elementMaps.psi_trace,
            self.u[0].femSpace.elementMaps.grad_psi_trace,
            self.elementBoundaryQuadratureWeights[('u', 0)],
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.psi_trace,
            self.u[0].femSpace.grad_psi_trace,
            self.u[0].femSpace.elementMaps.boundaryNormals,
            self.u[0].femSpace.elementMaps.boundaryJacobians,
            self.mesh.nElements_global,
            self.u[0].femSpace.dofMap.l2g,
            self.l2g[0]['freeGlobal'],
            self.mesh.elementDiametersArray,
            self.u[0].dof,
            self.coefficients.q_v,
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.coefficients.ebqe_v,
            self.numericalFlux.isDOFBoundary[0],
            self.numericalFlux.ebqe[('u', 0)],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.csrColumnOffsets_eb[(0, 0)],
            len(self.rowptr)-1,
            self.blendingFunctionDOFs, #alpha
            self.aux[0].femSpace.psi,
            self.aux[0].femSpace.grad_psi,
            self.dLow,
            self.coefficients.epsilon,
            self.coefficients.DISCRETE_UPWINDING,
            self.GALERKIN_SOLUTION,
            self.coefficients.PROBLEM_TYPE)

        # Load the Dirichlet conditions directly into residual
        if self.forceStrongConditions:
            scaling = 1.0
            for dofN in list(self.dirichletConditionsForceDOF.DOFBoundaryConditionsDict.keys()):
                gi = dofN
                for i in range(self.rowptr[gi], self.rowptr[gi + 1]):
                    gj = self.colind[i]
                    if (gj == gi):
                        self.nzval[i] = scaling
                    else:
                        self.nzval[i] = 0.0
                        #for k in range(self.rowptr[gj], self.rowptr[gj + 1]):
                        #    if (self.colind[k] == gi):
                        #        self.nzval[k] = 0.
        #

        logEvent("Jacobian ", level=10, data=jacobian)
        self.nonlinear_function_jacobian_evaluations += 1
        jacobian.fwrite("matdebug_p%s.vof.txt" % self.comm.rank())
        return jacobian

    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints, self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisHessianValuesRef(self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(self.timeIntegration.t, self.q)

        #mql. compute shape function of aux FE space
        self.aux[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.aux[0].femSpace.getBasisGradientValuesRef(self.elementQuadraturePoints)

    def calculateElementBoundaryQuadrature(self):
        pass

    def calculateExteriorElementBoundaryQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points on global element boundaries.

        This function should be called only when the mesh changes.
        """
        #
        # get physical locations of element boundary quadrature points
        #
        # assume all components live on the same mesh
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(self.elementBoundaryQuadraturePoints)
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
        pass

    def updateAfterMeshMotion(self):
        pass
