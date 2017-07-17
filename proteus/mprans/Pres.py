import proteus
import numpy
from proteus import *
from proteus.Transport import *
from proteus.Transport import OneLevelTransport
import os
from proteus import cfemIntegrals, Quadrature, Norms, Comm
from proteus.NonlinearSolvers import NonlinearEquation
from proteus.FemTools import (DOFBoundaryConditions,
                              FluxBoundaryConditions,
                              C0_AffineLinearOnSimplexWithNodalBasis)
from proteus.flcbdfWrappers import globalMax
from proteus.Profiling import memory
from proteus.Profiling import logEvent as log
from proteus.Transport import OneLevelTransport
from proteus.TransportCoefficients import TC_base
from proteus.SubgridError import SGE_base
from proteus.ShockCapturing import ShockCapturing_base
import cPres

class NumericalFlux(proteus.NumericalFlux.ConstantAdvection_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.ConstantAdvection_exterior.__init__(self,
                                                                  vt,getPointwiseBoundaryConditions,
                                                                  getAdvectiveFluxBoundaryConditions,
                                                                  getDiffusiveFluxBoundaryConditions)
        
class Coefficients(TC_base):
    r"""
    The coefficients for pressure solution

    Update is given by

    .. math::

       p^{k+1} - p^{k} - phi^{k+1} + \nabla\cdot(\mu \mathbf{u}^{k+1}) = 0
    """
    def __init__(self,
                 modelIndex=None,
                 fluidModelIndex=None,
                 pressureIncrementModelIndex=None,
                 useRotationalForm=False):
        """Construct a coefficients object

        :param pressureIncrementModelIndex: The index into the model list
        """
        TC_base.__init__(self,
                         nc = 1,
                         variableNames = ['p'],
                         reaction = {0:{0:'linear'}}, #  = p - p_last - phi
                         advection = {0:{0:'constant'}})# div  (\mu velocity)
        self.modelIndex = modelIndex
        self.fluidModelIndex = fluidModelIndex
        self.pressureIncrementModelIndex = pressureIncrementModelIndex
        self.useRotationalForm = useRotationalForm
    def attachModels(self,modelList):
        self.model = modelList[self.modelIndex]
        self.model.numericalFlux.ebqe[('u_last',0)] = self.model.ebqe[('u',0)].copy()
        self.model.q_p_sharp = self.model.q[('u',0)].copy()
        self.model.ebqe_p_sharp = self.model.ebqe[('u',0)].copy()
        self.model.q_grad_p_sharp = self.model.q[('grad(u)',0)].copy()
        self.model.ebqe_grad_p_sharp = self.model.ebqe[('grad(u)',0)].copy()
        self.model.u[0].femSpace.elementMaps.getBasisValuesRef(
            self.model.elementQuadraturePoints)
        self.model.u[0].femSpace.elementMaps.getBasisGradientValuesRef(
            self.model.elementQuadraturePoints)
        self.model.u[0].femSpace.getBasisValuesRef(self.model.elementQuadraturePoints)
        self.model.u[0].femSpace.getBasisGradientValuesRef(
            self.model.elementQuadraturePoints)
        self.model.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.model.elementBoundaryQuadraturePoints)
        self.model.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.model.elementBoundaryQuadraturePoints)
        self.model.u[0].femSpace.getBasisValuesTraceRef(
            self.model.elementBoundaryQuadraturePoints)
        self.model.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.model.elementBoundaryQuadraturePoints)
        self.pressureIncrementModel = modelList[
            self.pressureIncrementModelIndex]
        self.fluidModel = modelList[
            self.fluidModelIndex]
    def initializeMesh(self,mesh):
        """
        Give the TC object access to the mesh for any mesh-dependent information.
        """
        pass
    def initializeElementQuadrature(self,t,cq):
        """
        Give the TC object access to the element quadrature storage
        """
        cq[('u_last',0)] = cq[('u',0)].copy()
        self.q_massFlux = cq[('grad(u)',0)].copy()
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        """
        Give the TC object access to the element boundary quadrature storage
        """
        cebq[('u_last',0)] = cebq[('u',0)].copy()
        cebq_global[('u_last',0)] = cebq_global[('u',0)].copy()
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        """
        Give the TC object access to the exterior element boundary quadrature storage
        """
        cebqe[('u_last',0)] = cebqe[('u',0)].copy()
        self.ebqe_massFlux = cebqe[('grad(u)',0)].copy()
    def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
        """
        Give the TC object access to the generalized interpolation point storage. These points are used  to project nonlinear potentials (phi).
        """
        cip[('u_last',0)] = cip[('u',0)].copy()
    def preStep(self,t,firstStep=False):
        """
        Give the TC object an opportunity to modify itself before the time step.
        """
        self.q_massFlux[:] = self.fluidModel.q[('velocity',0)]
        np.multiply(self.fluidModel.coefficients.q_rho[:,:,np.newaxis],
                    self.q_massFlux,
                    out=self.q_massFlux)
        np.multiply(self.fluidModel.coefficients.q_nu[:,:,np.newaxis],
                    self.q_massFlux,
                    out=self.q_massFlux)
        self.ebqe_massFlux[:] = self.fluidModel.ebqe[('velocity',0)]
        np.multiply(self.fluidModel.coefficients.ebqe_rho[:,:,np.newaxis],
                    self.ebqe_massFlux,
                    out=self.ebqe_massFlux)
        np.multiply(self.fluidModel.coefficients.ebqe_nu[:,:,np.newaxis],
                    self.ebqe_massFlux,
                    out=self.ebqe_massFlux)
        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        """
        Give the TC object an opportunity to modify itself before the time step.
        """
        self.model.q[('u_last',0)][:] = self.model.q[('u',0)]
        self.model.ebqe[('u_last',0)][:] = self.model.ebqe[('u',0)]
        self.model.q_p_sharp[:] = self.model.q[('u',0)] + self.pressureIncrementModel.q[('u',0)]
        self.model.ebqe_p_sharp[:] = self.model.ebqe[('u',0)] + self.pressureIncrementModel.ebqe[('u',0)]
        self.model.q_grad_p_sharp[:] = self.model.q[('grad(u)',0)] + self.pressureIncrementModel.q[('grad(u)',0)]
        self.model.ebqe_grad_p_sharp[:] = self.model.ebqe[('grad(u)',0)] + self.pressureIncrementModel.ebqe[('grad(u)',0)]
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        self.evaluatePressure(t,c)
    def evaluatePressure(self,t,c):
        """
        Evaluate the coefficients after getting the specified velocity and density
        """
        # precompute the shapes to extract things we need from self.c_name[] dictionaries
        u_shape = c[('u',0)].shape
        grad_shape = c[('grad(u)',0)].shape
        if u_shape == self.pressureIncrementModel.q[('u',0)].shape:
            phi = self.pressureIncrementModel.q[('u',0)]
            rho = self.fluidModel.coefficients.q_rho
            nu = self.fluidModel.coefficients.q_nu
            velocity = self.fluidModel.q[('velocity',0)]
        elif u_shape == self.pressureIncrementModel.ebqe[('u',0)].shape:
            phi = self.pressureIncrementModel.ebqe[('u',0)]
            rho = self.fluidModel.coefficients.ebqe_rho
            nu = self.fluidModel.coefficients.ebqe_nu
            velocity = self.fluidModel.ebqe[('velocity',0)]
        # current and previous pressure values
        p = c[('u',0)]
        p_last = c[('u_last',0)]

        # set coefficients   p - p_star - phi + div (mu u) = 0
        #G&S11,p941,remark 5.5
        if self.useRotationalForm:
            for i in range(c[('f',0)].shape[-1]):
                c[('f',0)][...,i] = rho*nu*velocity[...,i]
        #G&S11,p92, eq 3.10
        c[('r',0)][:] = p - p_last - phi
        c[('dr',0,0)][:] = 1.0


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
                 movingDomain=False):
        from proteus import Comm
        #
        # set the objects describing the method and boundary conditions
        #
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
        self.ua = {}  # analytical solutions
        self.phi = phiDict
        self.dphi = {}
        self.matType = matType
        # mwf try to reuse test and trial information across components if
        # spaces are the same
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
        self.coefficients.initializeMesh(self.mesh)
        self.nc = self.coefficients.nc
        self.stabilization = stabilization
        self.shockCapturing = shockCapturing
        # no velocity post-processing for now
        self.conservativeFlux = conservativeFluxDict
        self.fluxBoundaryConditions = fluxBoundaryConditionsDict
        self.advectiveFluxBoundaryConditionsSetterDict = advectiveFluxBoundaryConditionsSetterDict
        self.diffusiveFluxBoundaryConditionsSetterDictDict = diffusiveFluxBoundaryConditionsSetterDictDict
        # determine whether  the stabilization term is nonlinear
        self.stabilizationIsNonlinear = False
        # cek come back
        if self.stabilization is not None:
            for ci in range(self.nc):
                if coefficients.mass.has_key(ci):
                    for flag in coefficients.mass[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.advection.has_key(ci):
                    for flag in coefficients.advection[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.diffusion.has_key(ci):
                    for diffusionDict in coefficients.diffusion[ci].values():
                        for flag in diffusionDict.values():
                            if flag != 'constant':
                                self.stabilizationIsNonlinear = True
                if coefficients.potential.has_key(ci):
                    for flag in coefficients.potential[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.reaction.has_key(ci):
                    for flag in coefficients.reaction[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
                if coefficients.hamiltonian.has_key(ci):
                    for flag in coefficients.hamiltonian[ci].values():
                        if flag == 'nonlinear':
                            self.stabilizationIsNonlinear = True
        # determine if we need element boundary storage
        self.elementBoundaryIntegrals = {}
        for ci in range(self.nc):
            self.elementBoundaryIntegrals[ci] = (
                (self.conservativeFlux is not None) or (
                    numericalFluxType is not None) or (
                    self.fluxBoundaryConditions[ci] == 'outFlow') or (
                    self.fluxBoundaryConditions[ci] == 'mixedFlow') or (
                    self.fluxBoundaryConditions[ci] == 'setFlow'))
        #
        # calculate some dimensions
        #
        # assume same space dim for all variables
        self.nSpace_global = self.u[0].femSpace.nSpace_global
        self.nDOF_trial_element = [
            u_j.femSpace.max_nDOF_element for u_j in self.u.values()]
        self.nDOF_phi_trial_element = [
            phi_k.femSpace.max_nDOF_element for phi_k in self.phi.values()]
        self.n_phi_ip_element = [
            phi_k.femSpace.referenceFiniteElement.interpolationConditions.nQuadraturePoints for phi_k in self.phi.values()]
        self.nDOF_test_element = [
            femSpace.max_nDOF_element for femSpace in self.testSpace.values()]
        self.nFreeDOF_global = [
            dc.nFreeDOF_global for dc in self.dirichletConditions.values()]
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
                if elementQuadrature.has_key(I):
                    elementQuadratureDict[I] = elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = elementQuadrature['default']
        else:
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = elementQuadrature
        if self.stabilization is not None:
            for I in self.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(I):
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature[I]
                    else:
                        elementQuadratureDict[
                            ('stab',) + I[1:]] = elementQuadrature['default']
                else:
                    elementQuadratureDict[
                        ('stab',) + I[1:]] = elementQuadrature
        if self.shockCapturing is not None:
            for ci in self.shockCapturing.components:
                if elemQuadIsDict:
                    if elementQuadrature.has_key(('numDiff', ci, ci)):
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            ('numDiff', ci, ci)]
                    else:
                        elementQuadratureDict[('numDiff', ci, ci)] = elementQuadrature[
                            'default']
                else:
                    elementQuadratureDict[
                        ('numDiff', ci, ci)] = elementQuadrature
        if massLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('m', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        if reactionLumping:
            for ci in self.coefficients.mass.keys():
                elementQuadratureDict[('r', ci)] = Quadrature.SimplexLobattoQuadrature(
                    self.nSpace_global, 1)
            for I in self.coefficients.elementIntegralKeys:
                elementQuadratureDict[
                    ('stab',) + I[1:]] = Quadrature.SimplexLobattoQuadrature(self.nSpace_global, 1)
        elementBoundaryQuadratureDict = {}
        if isinstance(elementBoundaryQuadrature, dict):  # set terms manually
            for I in self.coefficients.elementBoundaryIntegralKeys:
                if elementBoundaryQuadrature.has_key(I):
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[
                        I] = elementBoundaryQuadrature['default']
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
        (self.elementBoundaryQuadraturePoints, self.elementBoundaryQuadratureWeights,
         self.elementBoundaryQuadratureRuleIndeces) = Quadrature.buildUnion(elementBoundaryQuadratureDict)
        self.nElementBoundaryQuadraturePoints_elementBoundary = self.elementBoundaryQuadraturePoints.shape[
            0]
        self.nElementBoundaryQuadraturePoints_global = (
            self.mesh.nElements_global *
            self.mesh.nElementBoundaries_element *
            self.nElementBoundaryQuadraturePoints_elementBoundary)
        if type(self.u[0].femSpace) == C0_AffineLinearOnSimplexWithNodalBasis:
            # print self.nQuadraturePoints_element
            if self.nSpace_global == 3:
                assert(self.nQuadraturePoints_element == 5)
            elif self.nSpace_global == 2:
                assert(self.nQuadraturePoints_element == 6)
            elif self.nSpace_global == 1:
                assert(self.nQuadraturePoints_element == 3)

            # print self.nElementBoundaryQuadraturePoints_elementBoundary
            if self.nSpace_global == 3:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 2:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 4)
            elif self.nSpace_global == 1:
                assert(self.nElementBoundaryQuadraturePoints_elementBoundary == 1)
        #
        # simplified allocations for test==trial and also check if space is mixed or not
        #
        self.q = {}
        self.ebq = {}
        self.ebq_global = {}
        self.ebqe = {}
        self.phi_ip = {}
        # mesh
        #self.q['x'] = numpy.zeros((self.mesh.nElements_global,self.nQuadraturePoints_element,3),'d')
        self.ebqe['x'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             3),
            'd')
        self.q[('u', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')
        self.q[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global),
            'd')
        self.q[('r', 0)] = numpy.zeros(
            (self.mesh.nElements_global, self.nQuadraturePoints_element), 'd')

        self.ebqe[
            ('u',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[
            ('grad(u)',
             0)] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')

        self.points_elementBoundaryQuadrature = set()
        self.scalars_elementBoundaryQuadrature = set(
            [('u', ci) for ci in range(self.nc)])
        self.vectors_elementBoundaryQuadrature = set()
        self.tensors_elementBoundaryQuadrature = set()
        log(memory("element and element boundary Jacobians",
                   "OneLevelTransport"), level=4)
        self.inflowBoundaryBC = {}
        self.inflowBoundaryBC_values = {}
        self.inflowFlux = {}
        for cj in range(self.nc):
            self.inflowBoundaryBC[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,), 'i')
            self.inflowBoundaryBC_values[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global, self.nDOF_trial_element[cj]), 'd')
            self.inflowFlux[cj] = numpy.zeros(
                (self.mesh.nExteriorElementBoundaries_global,
                 self.nElementBoundaryQuadraturePoints_elementBoundary),
                'd')
        self.internalNodes = set(range(self.mesh.nNodes_global))
        # identify the internal nodes this is ought to be in mesh
        # \todo move this to mesh
        for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
            ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
            eN_global = self.mesh.elementBoundaryElementsArray[ebN, 0]
            ebN_element = self.mesh.elementBoundaryLocalElementBoundariesArray[
                ebN, 0]
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
        log("Updating local to global mappings", 2)
        self.updateLocal2Global()
        log("Building time integration object", 2)
        log(memory("inflowBC, internalNodes,updateLocal2Global",
                   "OneLevelTransport"), level=4)
        # mwf for interpolating subgrid error for gradients etc
        if self.stabilization and self.stabilization.usesGradientStabilization:
            self.timeIntegration = TimeIntegrationClass(
                self, integrateInterpolationPoints=True)
        else:
            self.timeIntegration = TimeIntegrationClass(self)

        if options is not None:
            self.timeIntegration.setFromOptions(options)
        log(memory("TimeIntegration", "OneLevelTransport"), level=4)
        log("Calculating numerical quadrature formulas", 2)
        self.calculateQuadrature()
        self.setupFieldStrides()

        comm = Comm.get()
        self.comm = comm
        if comm.size() > 1:
            assert numericalFluxType is not None and numericalFluxType.useWeakDirichletConditions, "You must use a numerical flux to apply weak boundary conditions for parallel runs"

        log(memory("stride+offset", "OneLevelTransport"), level=4)
        if numericalFluxType is not None:
            if options is None or options.periodicDirichletConditions is None:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict)
            else:
                self.numericalFlux = numericalFluxType(
                    self,
                    dofBoundaryConditionsSetterDict,
                    advectiveFluxBoundaryConditionsSetterDict,
                    diffusiveFluxBoundaryConditionsSetterDictDict,
                    options.periodicDirichletConditions)
        else:
            self.numericalFlux = None
        # set penalty terms
        # cek todo move into numerical flux initialization
        if self.ebq_global.has_key('penalty'):
            for ebN in range(self.mesh.nElementBoundaries_global):
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebq_global['penalty'][ebN, k] = self.numericalFlux.penalty_constant / (
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power)
        # penalty term
        # cek move  to Numerical flux initialization
        if self.ebqe.has_key('penalty'):
            for ebNE in range(self.mesh.nExteriorElementBoundaries_global):
                ebN = self.mesh.exteriorElementBoundariesArray[ebNE]
                for k in range(
                        self.nElementBoundaryQuadraturePoints_elementBoundary):
                    self.ebqe['penalty'][ebNE, k] = self.numericalFlux.penalty_constant / \
                        self.mesh.elementBoundaryDiametersArray[ebN]**self.numericalFlux.penalty_power
        log(memory("numericalFlux", "OneLevelTransport"), level=4)
        self.elementEffectiveDiametersArray = self.mesh.elementInnerDiametersArray
        #strong Dirichlet
        self.dirichletConditionsForceDOF= {0:DOFBoundaryConditions(self.u[cj].femSpace,dofBoundaryConditionsSetterDict[cj],weakDirichletConditions=False)}
        # use post processing tools to get conservative fluxes, None by default
        from proteus import PostProcessingTools
        self.velocityPostProcessor = PostProcessingTools.VelocityPostProcessingChooser(
            self)
        log(memory("velocity postprocessor", "OneLevelTransport"), level=4)
        # helper for writing out data storage
        from proteus import Archiver
        self.elementQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.elementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.exteriorElementBoundaryQuadratureDictionaryWriter = Archiver.XdmfWriter()
        self.globalResidualDummy = None
        compKernelFlag = 0
        self.pres = cPres.Pres(
            self.nSpace_global,
            self.nQuadraturePoints_element,
            self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
            self .u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
            self.testSpace[0].referenceFiniteElement.localFunctionSpace.dim,
            self.nElementBoundaryQuadraturePoints_elementBoundary,
            compKernelFlag)

    def calculateCoefficients(self):
        pass

    def calculateElementResidual(self):
        if self.globalResidualDummy is not None:
            self.getResidual(self.u[0].dof, self.globalResidualDummy)

    def getResidual(self, u, r):
        import pdb
        import copy
        from proteus.flcbdfWrappers import globalSum
        """
        Calculate the element residuals and add in to the global residual
        """
        r.fill(0.0)
        #set strong Dirichlet conditions
        for dofN,g in self.dirichletConditionsForceDOF[0].DOFBoundaryConditionsDict.iteritems():
            u[self.offset[0]+self.stride[0]*dofN] = g(self.dirichletConditionsForceDOF[0].DOFBoundaryPointDict[dofN],self.timeIntegration.t)#load the BC valu        # Load the unknowns into the finite element dof
        self.setUnknowns(u)

        # no flux boundary conditions
        self.pres.calculateResidual(  # element
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
            self.u[0].dof,
            self.q[('u', 0)],
            self.q[('grad(u)', 0)],
            self.q[('u_last',0)],
            self.coefficients.pressureIncrementModel.q[('u',0)],
            self.coefficients.q_massFlux,
            self.coefficients.ebqe_massFlux,
            self.ebqe[('u', 0)],
            self.ebqe[('grad(u)', 0)],
            self.offset[0], self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray)
        for dofN,g in self.dirichletConditionsForceDOF[0].DOFBoundaryConditionsDict.iteritems():
            r[self.offset[0]+self.stride[0]*dofN] = self.u[0].dof[dofN] - g(self.dirichletConditionsForceDOF[0].DOFBoundaryPointDict[dofN],self.timeIntegration.t)
        log("Global residual", level=9, data=r)
        self.nonlinear_function_evaluations += 1

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian, jacobian)
        self.pres.calculateJacobian(  # element
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
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            jacobian)
        for dofN in self.dirichletConditionsForceDOF[0].DOFBoundaryConditionsDict.keys():
            global_dofN = self.offset[0]+self.stride[0]*dofN
            self.nzval[numpy.where(self.colind == global_dofN)] = 0.0 #column
            self.nzval[self.rowptr[global_dofN]:self.rowptr[global_dofN+1]] = 0.0 #row
            zeroRow=True
            for i in range(self.rowptr[global_dofN],self.rowptr[global_dofN+1]):#row
                if (self.colind[i] == global_dofN):
                    self.nzval[i] = 1.0
                    zeroRow = False
            if zeroRow:
                raise RuntimeError("Jacobian has a zero row because sparse matrix has no diagonal entry at row "+`global_dofN`+". You probably need add diagonal mass or reaction term")
        log("Jacobian ", level=10, data=jacobian)
        self.nonlinear_function_jacobian_evaluations += 1
        return jacobian

    def calculateElementQuadrature(self):
        """
        Calculate the physical location and weights of the quadrature rules
        and the shape information at the quadrature points.

        This function should be called only when the mesh changes.
        """
        # self.u[0].femSpace.elementMaps.getValues(self.elementQuadraturePoints,
        #                                          self.q['x'])
        self.u[0].femSpace.elementMaps.getBasisValuesIP(
            self.u[0].femSpace.referenceFiniteElement.interpolationConditions.quadraturePointArray)
        self.u[0].femSpace.elementMaps.getBasisValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisValuesRef(self.elementQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesRef(
            self.elementQuadraturePoints)
        self.coefficients.initializeElementQuadrature(
            self.timeIntegration.t, self.q)
        if self.stabilization is not None:
            self.stabilization.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)
            self.stabilization.initializeTimeIntegration(self.timeIntegration)
        if self.shockCapturing is not None:
            self.shockCapturing.initializeElementQuadrature(
                self.mesh, self.timeIntegration.t, self.q)

    def calculateElementBoundaryQuadrature(self):
        pass

    def calculateExteriorElementBoundaryQuadrature(self):
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(
            self.elementBoundaryQuadraturePoints, self.ebqe['x'])
        # self.fluxBoundaryConditionsObjectsDict = dict([(cj, FluxBoundaryConditions(self.mesh,
        #                                                                            self.nElementBoundaryQuadraturePoints_elementBoundary,
        #                                                                            self.ebqe[('x')],
        #                                                                            self.advectiveFluxBoundaryConditionsSetterDict[cj],
        #                                                                            self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
        #                                                for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(
            self.timeIntegration.t, self.ebqe)

    def estimate_mt(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def updateAfterMeshMotion(self):
        pass


