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

class RKEV(proteus.TimeIntegration.SSP):
    from proteus import TimeIntegration
    """
    Wrapper for SSPRK time integration using EV

    ... more to come ...
    """

    def __init__(self, transport, timeOrder=1, runCFL=0.1, integrateInterpolationPoints=False):
        BackwardEuler.__init__(self, transport, integrateInterpolationPoints=integrateInterpolationPoints)
        self.trasport = transport
        self.runCFL = runCFL
        self.dtLast = None
        self.dtRatioMax = 2.0
        self.isAdaptive = True
        # About the cfl
        assert hasattr(transport, 'edge_based_cfl'), "No edge based cfl defined"
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
        maxCFL = max(maxCFL, globalMax(self.edge_based_cfl.max()))
        
        self.dt = old_div(self.runCFL, maxCFL)


        if self.transport.coefficients.fixed_dt is not None:
	    self.dt = self.transport.coefficients.fixed_dt
     
        if self.dtLast is None:
            self.dtLast = self.dt
        if old_div(self.dt, self.dtLast) > self.dtRatioMax:
            self.dt = self.dtLast * self.dtRatioMax
        self.t = self.tLast + self.dt
        self.substeps = [self.t for i in range(self.nStages)]  # Ignoring dif. time step levels 

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
                self.transport.u_dof_old[:] = self.u_dof_lstage[0]
                # update velocity at t=tn+dt for the second stage
                if self.transport.hasVelocityFieldAsFunction:
                    time = self.tLast + self.dt
                    self.transport.updateVelocityFieldAsFunction(time)
                    self.transport.updateVelocityAtDOFs(time)
                #
                # update uInlet at t=tn+dt for the second stage
                if self.transport.hasInletFunction:
                    time = self.tLast + self.dt
                    self.transport.update_uInlet_at_quad_points(time)
                #
            elif self.lstage == 2:
                logEvent("Second stage of SSP33 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                    self.u_dof_lstage[ci] *= old_div(1., 4.)
                    self.u_dof_lstage[ci] += 3. / 4. * self.u_dof_last[ci]
                # update u_dof_old
                self.transport.u_dof_old[:] = self.u_dof_lstage[0]
                # update velocity at t=tn+dt/2.0 for the third (and last) stage
                if self.transport.hasVelocityFieldAsFunction:
                    time = self.tLast + self.dt/2.0
                    self.transport.updateVelocityFieldAsFunction(time)
                    self.transport.updateVelocityAtDOFs(time)
                #
                # update uInlet at t=tn+dt/2.0 for the third (and last) stage
                if self.transport.hasInletFunction:
                    time = self.tLast + self.dt/2.0
                    self.transport.update_uInlet_at_quad_points(time)
                #
            else:
                logEvent("Third stage of SSP33 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                    self.u_dof_lstage[ci][:] *= old_div(2.0, 3.0)
                    self.u_dof_lstage[ci][:] += 1.0 / 3.0 * self.u_dof_last[ci]
                    # update solution to u[0].dof
                    self.transport.u[ci].dof[:] = self.u_dof_lstage[ci]
                # update u_dof_old
                self.transport.u_dof_old[:] = self.u_dof_last[0]
        elif self.timeOrder == 2:
            if self.lstage == 1:
                logEvent("First stage of SSP22 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                # Update u_dof_old
                self.transport.u_dof_old[:] = self.u_dof_lstage[0]
                # update velocity at t=tn+dt for the second (and last) stage
                if self.transport.hasVelocityFieldAsFunction:
                    time = self.tLast + self.dt
                    self.transport.updateVelocityFieldAsFunction(time)
                    self.transport.updateVelocityAtDOFs(time)
                #
                # update uInlet at t=tn+dt for the second (and last) stage
                if self.transport.hasInletFunction:
                    time = self.tLast + self.dt
                    self.transport.update_uInlet_at_quad_points(time)
                #
            else:
                logEvent("Second stage of SSP22 method finished", level=4)
                for ci in range(self.nc):
                    self.u_dof_lstage[ci][:] = self.transport.u[ci].dof
                    self.u_dof_lstage[ci][:] *= old_div(1., 2.)
                    self.u_dof_lstage[ci][:] += 1. / 2. * self.u_dof_last[ci]
                    # update solution to u[0].dof
                    self.transport.u[ci].dof[:] = self.u_dof_lstage[ci]
                # Update u_dof_old
                self.transport.u_dof_old[:] = self.u_dof_last[0]
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
    from proteus.ctransportCoefficients import VOFCoefficientsEvaluate
    from proteus.ctransportCoefficients import VolumeAveragedVOFCoefficientsEvaluate
    from proteus.cfemIntegrals import copyExteriorElementBoundaryValuesFromElementBoundaryValues

    def __init__(self,
                 GET_POINT_VALUES=1,
                 ME_model=0,
                 epsFact=0.0,
                 forceStrongConditions=0,
                 # OUTPUT quantDOFs
                 outputQuantDOFs=False,
                 periodicBCs=0, #Only works in 1D
                 rightBoundary=None, #for the periodicBCs
                 PROBLEM_TYPE=0,
                 ONE_DIM_PROBLEM=0,
                 PROJECT_INIT_CONDITION=0,
                 METHOD=4,
                 cE=1.0,
                 fixed_dt=None):

        self.cE=cE
        self.periodicBCs=periodicBCs
        self.rightBoundary=rightBoundary
        if self.periodicBCs==1:
            assert rightBoundary is not None, "provide a rightBoundary to use periodicBCs"
        self.fixed_dt=fixed_dt
        self.PROJECT_INIT_CONDITION=PROJECT_INIT_CONDITION
        self.METHOD=METHOD
        self.ONE_DIM_PROBLEM=ONE_DIM_PROBLEM
        # METHOD 4
        #   0: low-order
        #   1: high-order non-limited
        #   2: limiting without gamma indicator
        #   3: limiting with gamma indicator based on DK and CL
        #   4: limiting with gamma indicator based on MQL, DK and CK
        #
        self.GET_POINT_VALUES=GET_POINT_VALUES
        self.PROBLEM_TYPE=PROBLEM_TYPE
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

        time = self.model.timeIntegration.tLast
        if self.model.hasInletFunction:
            self.model.update_uInlet_at_quad_points(time)
        #
        
        # COMPUTE NEW VELOCITY (if given by user) #        
        if self.model.hasVelocityFieldAsFunction:
            self.model.updateVelocityFieldAsFunction(time)
            self.model.updateVelocityAtDOFs(time)
        #
        
        copyInstructions = {}
        return copyInstructions

    def postStep(self, t, firstStep=False):
        print "********************... ", self.model.u[0].dof.min(), self.model.u[0].dof.max()
        print "********************... ", len(self.model.u[0].dof)
        
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

        self.edge_based_cfl = numpy.zeros(self.u[0].dof.shape)
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

        # uInlet
        self.hasInletFunction = False
        if ('uInletFunction') in dir (options):
            self.uInletFunction = options.uInletFunction
            self.hasInletFunction = True
        #

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
        # For smoothness indicators #
        self.is_dof_external = numpy.zeros(self.u[0].dof.shape, 'd')
        self.is_dof_internal = numpy.zeros(self.u[0].dof.shape, 'd')
        self.is_boundary = numpy.zeros(self.u[0].dof.shape, 'd')
        self.is_cell_boundary = numpy.zeros(self.u[0].dof.shape, 'd')
        self.num_hi = numpy.zeros(self.u[0].dof.shape, 'd')
        self.den_hi = numpy.zeros(self.u[0].dof.shape, 'd')
        self.hi = numpy.zeros(self.u[0].dof.shape, 'd')
        self.gamma_dof = numpy.zeros(self.u[0].dof.shape, 'd')
        self.gamma_elem = numpy.zeros(self.mesh.nElements_global,'d')
        self.first_adjacent_dof_to_middle_dof = numpy.zeros(self.u[0].dof.shape, 'i')
        self.second_adjacent_dof_to_middle_dof = numpy.zeros(self.u[0].dof.shape, 'i')
        
        # mql. For edge based stabilization
        self.dLow=None
        self.quantDOFs = numpy.zeros(self.u[0].dof.shape, 'd')
        self.quantDOFs2 = numpy.zeros(self.u[0].dof.shape, 'd')
        self.quantDOFs3 = numpy.zeros(self.u[0].dof.shape, 'd')
        self.quantDOFs4 = numpy.zeros(self.u[0].dof.shape, 'd')
        self.quantDOFs5 = numpy.zeros(self.u[0].dof.shape, 'd')
        self.boundaryValues = None
        self.isBoundary = None

        # mql. For blending functions
        self.blendingFunctionDOFs = numpy.zeros(self.u[0].dof.shape,'d')
        # get Q2mesh_dof. This has the location of all DOFs in the Q2 mesh
        self.Q2mesh_dof = numpy.zeros((self.u[0].dof.size,3),'d')
        for eN in range(self.mesh.nElements_global):
            for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                gi = self.offset[0]+self.stride[0]*self.u[0].femSpace.dofMap.l2g[eN,i]
                self.Q2mesh_dof[gi][0] = self.u[0].femSpace.interpolationPoints[eN,i,0]
                self.Q2mesh_dof[gi][1] = self.u[0].femSpace.interpolationPoints[eN,i,1]
                self.Q2mesh_dof[gi][2] = self.u[0].femSpace.interpolationPoints[eN,i,2]

        # mql. METRICS #
        self.hasInitialCondition = False
        self.uInitial = None
        if ('uInitial') in dir (options):
            self.hasInitialCondition = True
            self.uInitial = options.uInitial
        #
            
        self.hasExactSolution = False
        self.exactSolution = None
        self.exactGrad = None
        if ('exactSolution') in dir(options):
            self.hasExactSolution = True
            self.exactSolution = options.exactSolution
            self.exactGrad = options.exactGrad
        #
        self.global_L1 = 0.0
        self.global_L2 = 0.0
        self.global_LInf = 0.0
        self.global_H1 = 0.0
        self.metricsAtEOS = open(self.name+"_metricsAtEOS.csv","w")
        self.metricsAtEOS.write('global_L1'+"\t"+
                                'global_L2'+"\t"+
                                'global_Linf'+"\t"
                                'global_H1'+"\n")

        self.cterm_global = None
        # Stuff for highOrderLim implementation #
        self.QL_MC_a = None
        self.QL_ML = None  # lumped mass matrix with low order space
        self.QL_MC_global = None  # consistent mass matrix with low order space

        self.QH_MC_a = None
        self.QH_ML = None  # lumped mass matrix with high order space
        self.QH_MC_global = None  # consistent mass matrix with high order space

        # for flux Jacobian at DOFs
        self.periodicBoundaryMap = None
        self.x = None
        self.y = None
        self.z = None
        self.u_vel_dofs = None
        self.v_vel_dofs = None

        self.intBernMat = None

        self.element_flux_i = None
        self.vVector = None

        self.xGradRHS = None
        self.yGradRHS = None
        self.xqNorm = None
        self.yqNorm = None
        self.qNorm = None

        self.umaxG = None
        self.uminG = None
        self.EntVisc = None
        self.xFlux_dof_old = None
        self.yFlux_dof_old = None
        self.uHDot = None
        self.QL_NNZ = 0
        self.Q1_sparsity=None
        self.uInlet_at_quad_points=None
        self.q_uInitial=None
        self.INIT_CONDITION_PROJECTED=False

        self.elem_patch_array = None
    #

    def getMetricsAtEOS(self):
        """
        Calculate some metrics at EOS (End Of Simulation)
        """
        time = self.timeIntegration.tLast
        nElements = self.q[('x')][:,:,0].shape[0]
        nQuad = self.q[('x')][:,:,0].shape[1]
        
        u_exact = numpy.zeros(self.q[('u',0)].shape,'d')
        gradx_u_exact = numpy.zeros(self.q[('u',0)].shape,'d')
        grady_u_exact = numpy.zeros(self.q[('u',0)].shape,'d')
        
        X = {0:self.q[('x')][:,:,0].ravel(),
             1:self.q[('x')][:,:,1].ravel(),
             2:self.q[('x')][:,:,2].ravel()}

        u_exact[:] = np.reshape(self.exactSolution[0](X,time),(nElements,nQuad))        
        gradx_u_exact[:] = np.reshape(self.exactGrad[0](X,time),(nElements,nQuad))
        grady_u_exact[:] = np.reshape(self.exactGrad[1](X,time),(nElements,nQuad))

        (global_L1,
         global_L2,
         global_LInf,
         global_H1) = self.blendedSpaces.calculateMetricsAtEOS(#element
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
             self.u[0].dof,
             u_exact,
             gradx_u_exact,
             grady_u_exact,
             self.offset[0],self.stride[0])

        from proteus.flcbdfWrappers import globalSum
        from proteus.flcbdfWrappers import globalMax
        # Interface metrics
        relative_norms=False
        self.global_L1 = globalSum(global_L1)
        self.global_L2 = globalSum(global_L2)
        self.global_LInf = globalMax(global_LInf)
        self.global_H1 = globalSum(global_H1)

        #####
        file_array=np.zeros((self.mesh.nElements_global,3),'d')
        file_array[:,0] = self.xElemCoord
        file_array[:,1] = self.yElemCoord
        file_array[:,2] = self.gamma_elem
        np.savetxt('si.txt',file_array)
        ####
    #
    
    def runAtEOS(self):
        # COMPUTE NORMS #
        if self.hasExactSolution:
            self.getMetricsAtEOS()
            self.metricsAtEOS.write(repr(self.global_L1)+"\t"+
                                    repr(np.sqrt(self.global_L2))+"\t"+
                                    repr(self.global_LInf)+"\t"+
                                    repr(np.sqrt(self.global_H1))+"\n")
            self.metricsAtEOS.flush()
        #
    #

    def calculateCoefficients(self):
        pass

    def getPeriodicBCs(self,rightBoundary):
	# NOTE: it only works in 1D
        arrayRightBoundary=[]
        self.periodicBoundaryMap=[]

        # find the indeces corresponding to the right boundary
        for i in range(len(self.x)):
            if self.x[i] == rightBoundary:
                arrayRightBoundary.append(i)
            #
        #
        for i in range(len(self.x)):        
            if self.x[i] == 0:
                yLeft = self.y[i]
                #print (i, self.x[i], self.y[i])
                for j in range(len(arrayRightBoundary)):
                    yRightIndex = arrayRightBoundary[j]
                    yRight = self.y[yRightIndex]
                    if yRight == yLeft:
                        self.periodicBoundaryMap.append([i,yRightIndex])
                        break
                    #
                #
            #
        #
    #
    
    def updateVelocityAtDOFs(self,time):
        if self.x is None:
            self.x = numpy.zeros(self.u[0].dof.shape,'d')
            self.y = numpy.zeros(self.u[0].dof.shape,'d')
            self.z = numpy.zeros(self.u[0].dof.shape,'d')
            for eN in range(self.mesh.nElements_global):
                for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                    gi = self.offset[0]+self.stride[0]*self.u[0].femSpace.dofMap.l2g[eN,i]
                    self.x[gi] = self.u[0].femSpace.interpolationPoints[eN,i,0]
                    self.y[gi] = self.u[0].femSpace.interpolationPoints[eN,i,1]
                    self.z[gi] = self.u[0].femSpace.interpolationPoints[eN,i,2]
                #
            #
        #
        X = {0:self.x,
             1:self.y,
             2:self.z}

        self.u_vel_dofs = np.zeros(len(self.x),'d')
        self.v_vel_dofs = np.zeros(len(self.x),'d')
        self.u_vel_dofs[:] = self.velocityFieldAsFunction[0](X,time)
        self.v_vel_dofs[:] = self.velocityFieldAsFunction[1](X,time)
    #

    def getIntBernMat(self):
        degree = self.u[0].femSpace.order if hasattr(self.u[0].femSpace,'order') else 1.0
        if degree==2:
            localCoord = {0: [-1.0,-1.0],
                          1: [-1.0, 1.0],
                          2: [ 1.0, 1.0],
                          3: [ 1.0,-1.0],
                          4: [-1.0, 0.0],
                          5: [ 0.0, 1.0],
                          6: [ 1.0, 0.0],
                          7: [ 0.0,-1.0],
                          8: [ 0.0, 0.0]}
        elif degree==3:
            localCoord = {0: [-1.0,-1.0],
                          1: [ 1.0,-1.0],
                          2: [ 1.0, 1.0],
                          3: [-1.0, 1.0],
                          #
                          4: [-1./3.,-1.0],
                          5: [ 1./3.,-1.0],
                          #
                          6: [ 1.0,-1./3.],
                          7: [ 1.0, 1./3.],
                          #
                          8: [ 1./3., 1.0],
                          9: [-1./3., 1.0],
                          #
                          10: [-1.0, 1./3.],
                          11: [-1.0,-1./3.],
                          #
                          12: [-1./3.,-1./3.],
                          13: [ 1./3.,-1./3.],
                          14: [-1./3., 1./3.],
                          15: [ 1./3., 1./3.]}
        #
        else:
            localCoord = {
                0: [-1.0,-1.0],
                1: [ 1.0,-1.0],
                2: [ 1.0, 1.0],
                3: [-1.0, 1.0],
                #
                4: [-1./2.,-1.0],
                5: [    0.,-1.0],
                6: [ 1./2.,-1.0],
                #
                7: [ 1.0,-1./2.],
                8: [ 1.0,    0.],
                9: [ 1.0, 1./2.],
                #
                10: [ 1./2., 1.0],
                11: [    0., 1.0],
                12: [-1./2., 1.0],
                #
                13: [-1.0, 1./2.],
                14: [-1.0,    0.],
                15: [-1.0,-1./2.],
                #
                16: [-1./2.,-1./2.],
                17: [    0.,-1./2.],
                18: [ 1./2.,-1./2.],
                #
                19: [-1./2., 0.],
                20: [    0., 0.],
                21: [ 1./2., 0.],
                #
                22: [-1./2., 1./2.],
                23: [    0., 1./2.],
                24: [ 1./2., 1./2.]}

        self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.basis[0](localCoord[0])
        self.intBernMat = np.zeros(len(self.Cx),'d')
        nDOF_test_element = len(self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim)
        for eN in range(self.mesh.nElements_global):
            for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                eN_i = eN*nDOF_test_element+i;
                #print ("***********... ",i)
                for j in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                    eN_i_j = eN_i*nDOF_test_element+j;
                    ij = self.csrRowIndeces[(0, 0)][eN,i] + self.csrColumnOffsets[(0, 0)][eN,i,j]
                    coord = localCoord[i]
                    v=self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.basis[j](coord)
                    self.intBernMat[ij] = v
                    #print (i,j, v, coord)
                #
            #
        #
        # check partition of unity
        #
        if False:
            for index in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                pu = 0
                coord = localCoord[index]
                for node in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                    pu+=self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.basis[node](coord)
                print (coord,pu)
            input("finished checking partition of unity")
        #
    #
    
    def update_uInlet_at_quad_points(self,time):
        # DO FINITE ELEMENT INTERPOLATION #
        logEvent("Computing boundary data at quad points", level=2)
        nElements = self.ebqe['x'][:,:,0].shape[0]
        nQuad = self.ebqe['x'][:,:,0].shape[1]
        
        ebqe_X = {0: self.ebqe['x'][:, :, 0].ravel(),
                  1: self.ebqe['x'][:, :, 1].ravel(),
                  2: self.ebqe['x'][:, :, 2].ravel()}
        
        self.uInlet_at_quad_points[:] = np.reshape(self.uInletFunction[0](ebqe_X,time),
                                                   (nElements,nQuad))
    #

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
            self.quantDOFs[:] = self.blendingFunctionDOFs
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
            self.quantDOFs[:] = self.blendingFunctionDOFs

    def updateVelocityFieldAsFunction(self,time):
        X = {0: self.q[('x')][:, :, 0],
             1: self.q[('x')][:, :, 1],
             2: self.q[('x')][:, :, 2]}
        self.coefficients.q_v[..., 0] = self.velocityFieldAsFunction[0](X, time)
        self.coefficients.q_v[..., 1] = self.velocityFieldAsFunction[1](X, time)
        if (self.nSpace_global == 3):
            self.coefficients.q_v[..., 2] = self.velocityFieldAsFunction[2](X, time)

        # BOUNDARY
        ebqe_X = {0: self.ebqe['x'][:, :, 0],
                  1: self.ebqe['x'][:, :, 1],
                  2: self.ebqe['x'][:, :, 2]}
        self.coefficients.ebqe_v[..., 0] = self.velocityFieldAsFunction[0](ebqe_X, time)
        self.coefficients.ebqe_v[..., 1] = self.velocityFieldAsFunction[1](ebqe_X, time)
        if (self.nSpace_global == 3):
            self.coefficients.ebqe_v[..., 2] = self.velocityFieldAsFunction[2](ebqe_X, time)

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
            #if x==0.0 or x==1.0 or y==0.0 or y==1.0:
            #    value=-1.0
            #    self.boundaryValues[gi] = value
            #    self.isBoundary[gi] = 1.0
            # INTERNAL BOUNDARY #
            #if (x >= 4./9-eps and x <= 5./9+eps and y >= 4./9-eps and y <= 5./9+eps):
            #    value=1.0
            #    self.boundaryValues[gi] = value
            #    self.isBoundary[gi] = 1.0
            if x==0:
                self.boundaryValues[gi] = -1.
                self.isBoundary[gi]=1.0
            elif x==1:
                self.boundaryValues[gi] = 1.
                self.isBoundary[gi]=1.0
        #
    #

    def compute_QH_lumped_mass_matrix(self):
        # assume a linear mass term
        dm = np.ones(self.q[('u', 0)].shape, 'd')
        QH_elementMassMatrix = np.zeros((self.mesh.nElements_global,
                                         self.nDOF_test_element[0],
                                         self.nDOF_trial_element[0]), 'd')
        # JACOBIANS (FOR ELEMENT TRANSFORMATION)
        q_J = np.zeros((self.mesh.nElements_global,
                        self.nQuadraturePoints_element,
                        self.nSpace_global,
                        self.nSpace_global),
                       'd')
        q_inverse_J = np.zeros((self.mesh.nElements_global,
                                self.nQuadraturePoints_element,
                                self.nSpace_global,
                                self.nSpace_global),
                               'd')
        q_det_J = np.zeros((self.mesh.nElements_global,
                            self.nQuadraturePoints_element),
                           'd')
        self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                         q_J,
                                                         q_inverse_J,
                                                         q_det_J)
        q_abs_det_J = np.abs(q_det_J)
        # SHAPE FUNCTIONS
        q_w = np.zeros((self.mesh.nElements_global,
                        self.nQuadraturePoints_element,
                        self.nDOF_test_element[0]),
                       'd')
        q_w_times_dV_m = q_w.copy()
        self.u[0].femSpace.getBasisValues(self.elementQuadraturePoints, q_w)
        cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u', 0)],
                                             q_abs_det_J,
                                             q_w,
                                             q_w_times_dV_m)
        cfemIntegrals.updateMassJacobian_weak_lowmem(dm,
                                                     q_w,
                                                     q_w_times_dV_m,
                                                     QH_elementMassMatrix)
        rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
        nnz = nzval.shape[-1]  # number of non-zero entries in sparse matrix
        self.QH_MC_a = nzval.copy()
        self.QH_MC_global = SparseMat(self.nFreeDOF_global[0],
                                      self.nFreeDOF_global[0],
                                      nnz,
                                      self.QH_MC_a,
                                      colind,
                                      rowptr)
        cfemIntegrals.zeroJacobian_CSR(self.nnz, self.QH_MC_global)
        cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                  self.l2g[0]['freeLocal'],
                                                                  self.l2g[0]['nFreeDOF'],
                                                                  self.l2g[0]['freeLocal'],
                                                                  self.csrRowIndeces[(0, 0)],
                                                                  self.csrColumnOffsets[(0, 0)],
                                                                  QH_elementMassMatrix,
                                                                  self.QH_MC_global)
        self.QH_ML = np.zeros((self.nFreeDOF_global[0],), 'd')
        for i in range(self.nFreeDOF_global[0]):
            self.QH_ML[i] = self.QH_MC_a[rowptr[i]:rowptr[i + 1]].sum()
        np.testing.assert_almost_equal(self.QH_ML.sum(),
                                       self.mesh.volume,
                                       err_msg="Trace of lumped mass matrix should be the domain volume", verbose=True)
        # END OF COMPUTING LUMPED MASS MATRIX #
        ##########################################################
        # COMPUTE (ELEMENT AND GLOBAL) PRECONDITIONED C MATRICES #
        ##########################################################
        QH_ML_eN = np.zeros(self.nDOF_test_element[0])
        self.prCxElem = np.zeros((self.mesh.nElements_global,
                                 self.nDOF_test_element[0],
                                 self.nDOF_trial_element[0]), 'd')
        self.prCyElem = np.zeros((self.mesh.nElements_global,
                                 self.nDOF_test_element[0],
                                 self.nDOF_trial_element[0]), 'd')
        self.prCTxElem = np.zeros((self.mesh.nElements_global,
                                  self.nDOF_test_element[0],
                                  self.nDOF_trial_element[0]), 'd')
        self.prCTyElem = np.zeros((self.mesh.nElements_global,
                                  self.nDOF_test_element[0],
                                  self.nDOF_trial_element[0]), 'd')
        for eN in range(len(QH_elementMassMatrix)):
            for i in range(self.nDOF_test_element[0]):
                QH_ML_eN[i] = sum(QH_elementMassMatrix[eN][i])
            #
            LeftPreconditioner_eN=np.dot(np.diag(QH_ML_eN),np.linalg.inv(QH_elementMassMatrix[eN]))
            RightPreconditioner_eN=np.dot(np.matrix.transpose(np.linalg.inv(QH_elementMassMatrix[eN])),np.diag(QH_ML_eN))
            self.prCxElem[eN] = np.dot(LeftPreconditioner_eN,self.cterm[0][eN])
            self.prCyElem[eN] = np.dot(LeftPreconditioner_eN,self.cterm[1][eN])
            #self.prCTxElem[eN]= self.prCxElem[eN].transpose()
            #self.prCTyElem[eN]= self.prCyElem[eN].transpose()
            self.prCTxElem[eN]= np.dot(self.cterm_transpose[0][eN],RightPreconditioner_eN)
            self.prCTyElem[eN]= np.dot(self.cterm_transpose[1][eN],RightPreconditioner_eN)
        #
        self.pr_cterm[0] = self.prCxElem
        self.pr_cterm[1] = self.prCyElem
        self.pr_cterm_transpose[0] = self.prCTxElem
        self.pr_cterm_transpose[1] = self.prCTyElem
        self.pr_cterm_a = {}
        self.pr_cterm_a_transpose = {}
        rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
        # assemble global C matrices
        for d in range(self.nSpace_global):  # spatial dimensions
            # assemble global preconditioned C matrices
            self.pr_cterm_a[d] = nzval.copy()
            self.pr_cterm_a_transpose[d] = nzval.copy()
            self.pr_cterm_global[d] = SparseMat(self.nFreeDOF_global[0],
                                                self.nFreeDOF_global[0],
                                                nnz,
                                                self.pr_cterm_a[d],
                                                colind,
                                                rowptr)
            self.pr_cterm_global_transpose[d] = SparseMat(self.nFreeDOF_global[0],
                                                          self.nFreeDOF_global[0],
                                                          nnz,
                                                          self.pr_cterm_a_transpose[d],
                                                          colind,
                                                          rowptr)
            cfemIntegrals.zeroJacobian_CSR(self.nnz, self.pr_cterm_global[d])
            cfemIntegrals.zeroJacobian_CSR(self.nnz, self.pr_cterm_global_transpose[d])
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0, 0)],
                                                                      self.csrColumnOffsets[(0, 0)],
                                                                      self.pr_cterm[d],
                                                                      self.pr_cterm_global[d])
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0, 0)],
                                                                      self.csrColumnOffsets[(0, 0)],
                                                                      self.pr_cterm_transpose[d],
                                                                      self.pr_cterm_global_transpose[d])
        #
        rowptr, colind, self.prCx = self.pr_cterm_global[0].getCSRrepresentation()
        rowptr, colind, self.prCy = self.pr_cterm_global[1].getCSRrepresentation()
        rowptr, colind, self.prCTx = self.pr_cterm_global_transpose[0].getCSRrepresentation()
        rowptr, colind, self.prCTy = self.pr_cterm_global_transpose[1].getCSRrepresentation()
        ###################################################################
        # END OF COMPUTING (ELEMENT AND GLOBAL) PRECONDITIONED C MATRICES #
        ###################################################################
    #
    
    def compute_QL_lumped_mass_matrix(self):
        # assume a linear mass term
        dm = np.ones(self.q[('u', 0)].shape, 'd')
        self.elementMassMatrix = np.zeros((self.mesh.nElements_global,
                                           self.nDOF_test_element[0],
                                           self.nDOF_trial_element[0]), 'd')
        # JACOBIANS (FOR ELEMENT TRANSFORMATION)
        q_J = np.zeros((self.mesh.nElements_global,
                        self.nQuadraturePoints_element,
                        self.nSpace_global,
                        self.nSpace_global),
                       'd')
        q_inverse_J = np.zeros((self.mesh.nElements_global,
                                self.nQuadraturePoints_element,
                                self.nSpace_global,
                                self.nSpace_global),
                               'd')
        q_det_J = np.zeros((self.mesh.nElements_global,
                            self.nQuadraturePoints_element),
                           'd')
        self.aux[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                           q_J,
                                                           q_inverse_J,
                                                           q_det_J)
        q_abs_det_J = np.abs(q_det_J)
        # SHAPE FUNCTIONS
        q_w = np.zeros((self.mesh.nElements_global,
                        self.nQuadraturePoints_element,
                        self.nDOF_test_element[0]),
                       'd')
        q_w_times_dV_m = q_w.copy()
        self.aux[0].femSpace.getBasisValues(self.elementQuadraturePoints, q_w)
        cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u', 0)],
                                             q_abs_det_J,
                                             q_w,
                                             q_w_times_dV_m)
        cfemIntegrals.updateMassJacobian_weak_lowmem(dm,
                                                     q_w,
                                                     q_w_times_dV_m,
                                                     self.elementMassMatrix)
        rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
        nnz = nzval.shape[-1]  # number of non-zero entries in sparse matrix
        self.QL_MC_a = nzval.copy()
        self.QL_MC_global = SparseMat(self.nFreeDOF_global[0],
                                      self.nFreeDOF_global[0],
                                      nnz,
                                      self.QL_MC_a,
                                      colind,
                                      rowptr)
        cfemIntegrals.zeroJacobian_CSR(self.nnz, self.QL_MC_global)
        cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                  self.l2g[0]['freeLocal'],
                                                                  self.l2g[0]['nFreeDOF'],
                                                                  self.l2g[0]['freeLocal'],
                                                                  self.csrRowIndeces[(0, 0)],
                                                                  self.csrColumnOffsets[(0, 0)],
                                                                  self.elementMassMatrix,
                                                                  self.QL_MC_global)
        self.QL_ML = np.zeros((self.nFreeDOF_global[0],), 'd')
        self.Q1_sparsity = np.zeros(self.Cx.shape,'d')
        ij=0
        for i in range(self.nFreeDOF_global[0]):
            self.QL_ML[i] = self.QL_MC_a[rowptr[i]:rowptr[i + 1]].sum()
            for j in range(self.rowptr[i], self.rowptr[i+1]):
                if self.QL_MC_a[ij]>1E-10:
                    self.Q1_sparsity[ij]=1.0
                    self.QL_NNZ+=1
                #
                ij+=1
            #
        #
        
        np.testing.assert_almost_equal(self.QL_ML.sum(),
                                       self.mesh.volume,
                                       err_msg="Trace of lumped mass matrix should be the domain volume", verbose=True)
        # END OF COMPUTING LUMPED MASS MATRIX #    
        ######################################
        # COMPUTE THE INVERSE OF ML_minus_MC # using the linear space
        ######################################
        element_ML_minus_MC = np.zeros((self.mesh.nElements_global,
                                        self.nDOF_test_element[0],
                                        self.nDOF_trial_element[0]), 'd')
        self.inv_element_ML_minus_MC = np.zeros((self.mesh.nElements_global,
                                                 self.nDOF_test_element[0],
                                                 self.nDOF_trial_element[0]), 'd')
        for eN in range(len(self.elementMassMatrix)):
            for i in range(self.nDOF_test_element[0]):
                lumpedElementMassMatrix = sum(self.elementMassMatrix[eN][i])
                for j in range(self.nDOF_trial_element[0]):
                    if i==j:
                        element_ML_minus_MC[eN][i][j] = lumpedElementMassMatrix - self.elementMassMatrix[eN][i][j]
                    else:
                        element_ML_minus_MC[eN][i][j] = - self.elementMassMatrix[eN][i][j]
                    #
                #
            #
            element_ML_minus_MC[eN][0][:]=0
            element_ML_minus_MC[eN][0][0]=1
            self.inv_element_ML_minus_MC[eN] = np.linalg.inv(element_ML_minus_MC[eN])
        ###############################################
        # END OF COMPUTING THE INVERSE OF ML_minus_MC #
        ###############################################
        
    def compute_c_matrices(self,use_Q1_lagrange=False):
        self.cterm = {}
        self.cterm_a = {}
        self.cterm_global = {}
        self.cterm_transpose = {}        
        self.cterm_a_transpose = {}
        self.cterm_global_transpose = {}
        #
        self.pr_cterm_global = {}
        self.pr_cterm_global_transpose = {}
        self.pr_cterm = {}
        self.pr_cterm_transpose = {}
        #
        
        rowptr, colind, nzval = self.jacobian.getCSRrepresentation()
        nnz = nzval.shape[-1]  # number of non-zero entries in sparse matrix
        di = self.q[('grad(u)', 0)].copy()  # direction of derivative
        # JACOBIANS (FOR ELEMENT TRANSFORMATION)
        q_J = np.zeros((self.mesh.nElements_global,
                        self.nQuadraturePoints_element,
                        self.nSpace_global,
                        self.nSpace_global),
                       'd')
        q_inverse_J = np.zeros((self.mesh.nElements_global,
                                self.nQuadraturePoints_element,
                                self.nSpace_global,
                                self.nSpace_global),
                               'd')
        q_det_J = np.zeros((self.mesh.nElements_global,
                            self.nQuadraturePoints_element),
                           'd')
        if use_Q1_lagrange:
            self.aux[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                               q_J,
                                                               q_inverse_J,
                                                               q_det_J)
        else:
            self.u[0].femSpace.elementMaps.getJacobianValues(self.elementQuadraturePoints,
                                                             q_J,
                                                             q_inverse_J,
                                                             q_det_J)
        #
        q_abs_det_J = np.abs(q_det_J)
        # SHAPE FUNCTIONS
        q_w = np.zeros((self.mesh.nElements_global,
                        self.nQuadraturePoints_element,
                        self.nDOF_test_element[0]),
                       'd')
        q_w_times_dV_m = q_w.copy()
        if use_Q1_lagrange:
            self.aux[0].femSpace.getBasisValues(self.elementQuadraturePoints, q_w)
        else:
            self.u[0].femSpace.getBasisValues(self.elementQuadraturePoints, q_w)
        cfemIntegrals.calculateWeightedShape(self.elementQuadratureWeights[('u', 0)],
                                             q_abs_det_J,
                                             q_w,
                                             q_w_times_dV_m)
        # GRADIENT OF TEST FUNCTIONS
        q_grad_w = np.zeros((self.mesh.nElements_global,
                             self.nQuadraturePoints_element,
                             self.nDOF_test_element[0],
                             self.nSpace_global),
                            'd')
        if use_Q1_lagrange:
            self.aux[0].femSpace.getBasisGradientValues(self.elementQuadraturePoints,
                                                        q_inverse_J,
                                                        q_grad_w)
        else:
            self.u[0].femSpace.getBasisGradientValues(self.elementQuadraturePoints,
                                                      q_inverse_J,
                                                      q_grad_w)
        q_grad_w_times_dV_f = np.zeros((self.mesh.nElements_global,
                                        self.nQuadraturePoints_element,
                                        self.nDOF_test_element[0],
                                        self.nSpace_global),
                                       'd')
        cfemIntegrals.calculateWeightedShapeGradients(self.elementQuadratureWeights[('u', 0)],
                                                      q_abs_det_J,
                                                      q_grad_w,
                                                      q_grad_w_times_dV_f)
        for d in range(self.nSpace_global):  # spatial dimensions
            # C matrices
            self.cterm[d] = np.zeros((self.mesh.nElements_global,
                                      self.nDOF_test_element[0],
                                      self.nDOF_trial_element[0]), 'd')
            self.cterm_a[d] = nzval.copy()
            self.cterm_global[d] = SparseMat(self.nFreeDOF_global[0],
                                             self.nFreeDOF_global[0],
                                             nnz,
                                             self.cterm_a[d],
                                             colind,
                                             rowptr)
            cfemIntegrals.zeroJacobian_CSR(self.nnz, self.cterm_global[d])
            di[:] = 0.0
            di[..., d] = 1.0
            cfemIntegrals.updateHamiltonianJacobian_weak_lowmem(di,
                                                                q_grad_w_times_dV_f,
                                                                q_w,
                                                                self.cterm[d])  # int[(di*grad(wj))*wi*dV]
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0, 0)],
                                                                      self.csrColumnOffsets[(0, 0)],
                                                                      self.cterm[d],
                                                                      self.cterm_global[d])
            # C Transpose matrices
            self.cterm_transpose[d] = np.zeros((self.mesh.nElements_global,
                                                self.nDOF_test_element[0],
                                                self.nDOF_trial_element[0]), 'd')
            self.cterm_a_transpose[d] = nzval.copy()
            self.cterm_global_transpose[d] = SparseMat(self.nFreeDOF_global[0],
                                                       self.nFreeDOF_global[0],
                                                       nnz,
                                                       self.cterm_a_transpose[d],
                                                       colind,
                                                       rowptr)
            cfemIntegrals.zeroJacobian_CSR(self.nnz, self.cterm_global_transpose[d])
            di[:] = 0.0
            di[..., d] = -1.0
            cfemIntegrals.updateAdvectionJacobian_weak_lowmem(di,
                                                              q_w,
                                                              q_grad_w_times_dV_f,
                                                              self.cterm_transpose[d])  # -int[(-di*grad(wi))*wj*dV]
            cfemIntegrals.updateGlobalJacobianFromElementJacobian_CSR(self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.l2g[0]['nFreeDOF'],
                                                                      self.l2g[0]['freeLocal'],
                                                                      self.csrRowIndeces[(0, 0)],
                                                                      self.csrColumnOffsets[(0, 0)],
                                                                      self.cterm_transpose[d],
                                                                      self.cterm_global_transpose[d])
        #
        rowptr, colind, self.Cx = self.cterm_global[0].getCSRrepresentation()
        rowptr, colind, self.Cy = self.cterm_global[1].getCSRrepresentation()
        rowptr, colind, self.CTx = self.cterm_global_transpose[0].getCSRrepresentation()
        rowptr, colind, self.CTy = self.cterm_global_transpose[1].getCSRrepresentation()
    #

    def compute_element_patch(self):
        DoFs_in_elem = np.empty((self.mesh.nElements_global,),dtype=object)
        elem_patch = np.empty((self.mesh.nElements_global,),dtype=object)
        for index in range(len(DoFs_in_elem)):
            DoFs_in_elem[index]=[]
            elem_patch[index]=[]
        #
        elem_in_DoFs = np.empty((len(self.u[0].dof),),dtype=object)
        for index in range(len(elem_in_DoFs)):
            elem_in_DoFs[index]=[]
        #

        # loop on elements
        for eN in range(self.mesh.nElements_global):
            for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                offset_u = self.offset[0]; stride_u = self.stride[0]
                l2g = self.l2g[0]['freeGlobal']
                gi = offset_u + stride_u*l2g[eN,i]
                DoFs_in_elem[eN].append(gi)
                elem_in_DoFs[gi].append(eN)
            #
        #
        
        for eN in range(self.mesh.nElements_global):
            elem_patch[eN].append(eN)
            for DoF in DoFs_in_elem[eN]:
                for elem in elem_in_DoFs[DoF]:
                    if (elem in elem_patch[eN]) == False:
                        elem_patch[eN].append(elem)
                #
            #
        #

        self.elem_patch_array = -np.ones((self.mesh.nElements_global,9),'i')
        for eN,elem_list in enumerate(elem_patch):
            for j,elem in enumerate(elem_list):
                self.elem_patch_array[eN,j] = elem
            #
        #

        # get coordinates of the center point of the elements
        self.xElemCoord = np.zeros(self.mesh.nElements_global,'d')
        self.yElemCoord = np.zeros(self.mesh.nElements_global,'d')
        if self.x is None:
            self.x = numpy.zeros(self.u[0].dof.shape,'d')
            self.y = numpy.zeros(self.u[0].dof.shape,'d')
            self.z = numpy.zeros(self.u[0].dof.shape,'d')
            for eN in range(self.mesh.nElements_global):
                for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                    gi = self.offset[0]+self.stride[0]*self.u[0].femSpace.dofMap.l2g[eN,i]
                    self.x[gi] = self.u[0].femSpace.interpolationPoints[eN,i,0]
                    self.y[gi] = self.u[0].femSpace.interpolationPoints[eN,i,1]
                    self.z[gi] = self.u[0].femSpace.interpolationPoints[eN,i,2]
                #
            #
        #
        for eN in range(self.mesh.nElements_global):
            for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                if i==8:
                    offset_u = self.offset[0]; stride_u = self.stride[0]
                    l2g = self.l2g[0]['freeGlobal']
                    gi = offset_u + stride_u*l2g[eN,i]
                    self.xElemCoord[eN] = self.x[gi]
                    self.yElemCoord[eN] = self.y[gi]
        
    #
        
    def compute_matrices(self):
        if self.cterm_global is None:
            self.dLow = np.zeros(self.nnz,'d')
            self.compute_c_matrices()
            self.compute_QL_lumped_mass_matrix()
            self.compute_QH_lumped_mass_matrix()
            if self.coefficients.GET_POINT_VALUES==1:
                self.getIntBernMat()
            #
            self.uInlet_at_quad_points=numpy.zeros((self.mesh.nExteriorElementBoundaries_global,
                                                    self.nElementBoundaryQuadraturePoints_elementBoundary), 'd')
            if self.hasInletFunction:
                self.update_uInlet_at_quad_points(time=0)
            #
            if self.hasVelocityFieldAsFunction:
                self.updateVelocityFieldAsFunction(time=0)
                self.updateVelocityAtDOFs(time=0)
            else:
                self.u_vel_dofs = np.zeros(self.u[0].dof.shape,'d')
                self.v_vel_dofs = np.zeros(self.u[0].dof.shape,'d')
            #
            self.element_flux_i = np.zeros((self.mesh.nElements_global,
                                            self.nDOF_test_element[0]),'d')
            self.vVector = np.zeros((self.mesh.nElements_global,
                                     self.nDOF_test_element[0]),'d')
            self.element_flux_qij = np.zeros((self.mesh.nElements_global,
                                              self.nDOF_test_element[0],
                                              self.nDOF_trial_element[0]), 'd')
            self.flux_qij = np.zeros(len(self.Cx),'d')
            self.dLowElem = np.zeros((self.mesh.nElements_global,
                                      self.nDOF_test_element[0],
                                      self.nDOF_trial_element[0]), 'd')
            self.xGradRHS = numpy.zeros(self.u[0].dof.shape, 'd')
            self.yGradRHS = numpy.zeros(self.u[0].dof.shape, 'd')
            self.xqNorm = numpy.zeros(self.u[0].dof.shape, 'd')
            self.yqNorm = numpy.zeros(self.u[0].dof.shape, 'd')
            self.qNorm = numpy.zeros(self.u[0].dof.shape, 'd')

            self.EntVisc = numpy.zeros(self.u[0].dof.shape, 'd')
            self.uHDot = numpy.zeros(self.u[0].dof.shape, 'd')

            self.xFlux_dof_old = numpy.zeros(self.u[0].dof.shape, 'd')
            self.yFlux_dof_old = numpy.zeros(self.u[0].dof.shape, 'd')
            
            self.calculateResidual = self.blendedSpaces.calculateResidual
            #self.calculateResidual = self.blendedSpaces.calculateResidualEntropyVisc
            self.calculateJacobian = self.blendedSpaces.calculateMassMatrix

            self.q_uInitial=numpy.zeros(self.q[('u',0)].shape,'d')
            self.getAdjacent_DoFs_to_middle_dof()
            
            ## FOR DEBUGGING ##
            # compare lumped mass matrices ? #
            if False:
                for i in range(self.nFreeDOF_global[0]):
                    print self.QL_ML[i], self.QH_ML[i]
                print "TOTAL MASS (low-order, high-order): ", self.QL_ML.sum(), self.QH_ML.sum()
            #
            # Print matrix? #
            if False:
                ij=0
                for i in range(self.nFreeDOF_global[0]):
                    for j in range(self.nFreeDOF_global[0]):
                        print "%8.4e" % self.Cx[ij], '\t',
                        ij+=1
                    #
                    print 
                #
            #
            # END OF DEBUGGING #
        #
    #

    def getAdjacent_DoFs_to_middle_dof(self):
        for eN in range(self.mesh.nElements_global):
            for i in self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.range_dim:
                offset_u = self.offset[0]; stride_u = self.stride[0]
                l2g = self.l2g[0]['freeGlobal']
                gi = offset_u + stride_u*l2g[eN,i]
                
                self.den_hi[gi] += 1
                
                if i<4:
                    self.is_dof_external[gi] = 1;
                elif i==8:
                    self.is_dof_internal[gi] = 1;
                else: # i>=4 and i<=7
		    gj1 = 0
                    gj2 = 0
		    if i==4:
			gj1 = offset_u+stride_u*l2g[eN,0] 
			gj2 = offset_u+stride_u*l2g[eN,1] 
		    elif i==5:
			gj1 = offset_u+stride_u*l2g[eN,1] 
			gj2 = offset_u+stride_u*l2g[eN,2] 
		    elif i==6:
			gj1 = offset_u+stride_u*l2g[eN,2] 
			gj2 = offset_u+stride_u*l2g[eN,3] 
		    else:
			gj1 = offset_u+stride_u*l2g[eN,3] 
			gj2 = offset_u+stride_u*l2g[eN,0] 
                    #
		    self.first_adjacent_dof_to_middle_dof[gi]  = gj1;
		    self.second_adjacent_dof_to_middle_dof[gi] = gj2;
                #
            #
        #
        
    def getResidual(self, u, r):
        import pdb
        import copy
        
        ## COMPUTE C MATRICES ##
        self.compute_matrices()

        if self.elem_patch_array is None:
            self.compute_element_patch()
        
        if self.coefficients.periodicBCs and self.periodicBoundaryMap is None:
            self.getPeriodicBCs(rightBoundary=self.coefficients.rightBoundary)
        """
        Calculate the element residuals and add in to the global residual
        """
        if self.u_dof_old is None:
            # Pass initial condition to u_dof_old
            self.u_dof_old = numpy.copy(self.u[0].dof)
            self.umaxG = self.u_dof_old.max()
            self.uminG = self.u_dof_old.min()
        #
        
        # Load the unknowns into the finite element dof
        self.timeIntegration.calculateCoefs()
        self.timeIntegration.calculateU(u)
        self.setUnknowns(self.timeIntegration.u)

        # init to zero some vectors
        #self.dLow.fill(0.0)
        self.vVector.fill(0.0)
        r.fill(0.0)
        self.flux_qij.fill(0.0)
        self.element_flux_i.fill(0.0)
        self.edge_based_cfl.fill(0.0)

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
            self.u[0].femSpace.Hessian_psi,
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
            self.aux[0].femSpace.psi_trace,
            self.dLow,
            self.coefficients.PROBLEM_TYPE,
            self.coefficients.ONE_DIM_PROBLEM,
            self.coefficients.METHOD,
            self.quantDOFs,
            self.quantDOFs4,
            self.quantDOFs5,
            # For highOrderLim
            self.q_uInitial,
            self.uInlet_at_quad_points,
            self.coefficients.GET_POINT_VALUES,
            self.flux_qij,
            self.element_flux_qij,
            self.elementMassMatrix,
            self.vVector,
            self.element_flux_i,
            self.intBernMat,
            self.edge_based_cfl,
            # Flux Jacobian at DOFs
            self.u_vel_dofs,
            self.v_vel_dofs,
            # lumped mass matrices
            #self.QL_ML,
            self.QH_ML,
            # inverse of dissipative mass matrix
            self.inv_element_ML_minus_MC,
            self.xFlux_dof_old,
            self.yFlux_dof_old,
            self.umaxG,
            self.uminG,
            self.uHDot,
            self.coefficients.cE,
            self.EntVisc,
            # for smoothness indicator
            self.elem_patch_array,
            self.gamma_elem,
            self.gamma_dof,
            self.first_adjacent_dof_to_middle_dof,
            self.second_adjacent_dof_to_middle_dof,
            self.is_dof_external,
	    self.is_dof_internal,
	    self.den_hi,
            # C-matrices
            self.Cx,
            self.Cy,
            self.CTx,
            self.CTy,
            self.prCx,
            self.prCy,
            self.prCTx,
            self.prCTy,
            self.cterm[0], #CxElem
            self.cterm[1], #CyElem
            self.cterm_transpose[0], #CTxElem
            self.cterm_transpose[1], #CTyElem
            self.pr_cterm[0], #prCxElem
            self.pr_cterm[1], #prCyElem
            self.pr_cterm_transpose[0], #prCTxElem
            self.pr_cterm_transpose[1], #prCTyElem
            self.dLowElem,
            self.Q1_sparsity,
            self.qNorm,
            self.xqNorm,
            self.yqNorm,
            self.xGradRHS,
            self.yGradRHS)

        #self.gamma_dof[:]=1.0
        self.quantDOFs2[:] = self.EntVisc
        #self.quantDOFs3[:] = self.xGradRHS
        #self.quantDOFs4[:] = self.yGradRHS
        self.quantDOFs3[:] = self.gamma_dof
        #self.quantDOFs4[:] = self.xGradRHS
        #self.quantDOFs5[:] = self.yGradRHS
        self.quantDOFs5[:] = self.qNorm
        
        #import pdb; pdb.set_trace()
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
            self.dLow)

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
        self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
                                                 self.q['x'])
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

        #mql. compute shape function (trace) of aux FE space
        self.aux[0].femSpace.getBasisValuesTraceRef(self.elementBoundaryQuadraturePoints)

    def estimate_mt(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass

    def updateAfterMeshMotion(self):
        pass
