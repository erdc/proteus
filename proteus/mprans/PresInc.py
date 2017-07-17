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
import cPresInc
class NumericalFlux(proteus.NumericalFlux.ConstantAdvection_Diffusion_SIPG_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.ConstantAdvection_Diffusion_SIPG_exterior.__init__(self,vt,getPointwiseBoundaryConditions,
                                                                                        getAdvectiveFluxBoundaryConditions,
                                                                                        getDiffusiveFluxBoundaryConditions)


class Coefficients(TC_base):
    r"""
    The coefficients for pressure increment solution

    Update is given by

    .. math::

       \nabla\cdot( -a \nabla \phi^{k+1} - \mathbf{q^t}^{k+1} ) = 0
       a = \frac{\tau (1-\theta_s)}{\rho_f} + \frac{\tau \theta_s}{\rho_s}
       q^t = (1-\theta_s) v_f + \theta_s v_s
    """
    def __init__(self,
                 rho_f_min=998.0,
                 rho_s_min=998.0,
                 nd=2,
                 modelIndex = None,
                 fluidModelIndex = None):
        """Construct a coefficients object

        :param modelIndex: This model's index into the model list
        :param fluidModelIndex: The fluid momentum model's index
        """
        assert(nd in [2,3])
        self.nd = nd
        if self.nd == 2:
            sdInfo    = {(0,0):(np.array([0,1,2],dtype='i'),
                                np.array([0,1],dtype='i'))}
        else:
            sdInfo    = {(0,0):(np.array([0,1,2,3],dtype='i'),
                                np.array([0,1,2],dtype='i'))}
        TC_base.__init__(self,
                         nc = 1,
                         variableNames = ['pInc'],
                         diffusion = {0:{0:{0:'constant'}}},
                         potential = {0:{0:'u'}},
                         advection = {0:{0:'constant'}},
                         sparseDiffusionTensors=sdInfo,
                         useSparseDiffusion = True)
        self.rho_f_min = rho_f_min
        self.rho_s_min = rho_s_min
        self.modelIndex = modelIndex
        self.fluidModelIndex = fluidModelIndex

    def attachModels(self,modelList):
        """
        Attach the model for velocity and density to PresureIncrement model
        """
        self.model = modelList[self.modelIndex]
        self.fluidModel = modelList[self.fluidModelIndex]

    def initializeMesh(self,mesh):
        """
        Give the TC object access to the mesh for any mesh-dependent information.
        """
        pass

    def initializeElementQuadrature(self,t,cq):
        """
        Give the TC object access to the element quadrature storage
        """
        pass
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        """
        Give the TC object access to the element boundary quadrature storage
        """
        pass
    def initializeGlobalExteriorElementBoundaryQuadrature(self,t,cebqe):
        """
        Give the TC object access to the exterior element boundary quadrature storage
        """
        pass
    def initializeGeneralizedInterpolationPointQuadrature(self,t,cip):
        """
        Give the TC object access to the generalized interpolation point storage. These points are used  to project nonlinear potentials (phi).
        """
        pass
    def preStep(self,t,firstStep=False):
        """
        Move the current values to values_last to keep cached set of values for bdf1 algorithm
        """
        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        """
        Update the fluid velocities
        """
        alphaBDF = self.fluidModel.timeIntegration.alpha_bdf
        for i in range(self.fluidModel.q[('velocity',0)].shape[-1]):
            self.fluidModel.q[('velocity',0)][...,i] -= self.model.q[('grad(u)',0)][...,i]/(self.rho_f_min*alphaBDF)
            #cek hack, need to do scale this right for 3p flow
            self.fluidModel.ebqe[('velocity',0)][...,i] = (self.model.ebqe[('advectiveFlux',0)]+self.model.ebqe[('diffusiveFlux',0,0)])*self.model.ebqe['n'][...,i]
            self.fluidModel.coefficients.q_velocity_solid[...,i] -= self.model.q[('grad(u)',0)][...,i]/(self.rho_s_min*alphaBDF)
            self.fluidModel.coefficients.ebqe_velocity_solid[...,i] -= self.model.ebqe[('grad(u)',0)][...,i]/(self.rho_s_min*alphaBDF)
        self.fluidModel.stabilization.v_last[:] = self.fluidModel.q[('velocity',0)]
        self.fluidModel.coefficients.ebqe_velocity_last[:] = self.fluidModel.ebqe[('velocity',0)]
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        self.evaluatePresureIncrement(t,c)
    def evaluatePresureIncrement(self,t,c):
        """
        Evaluate the coefficients after getting the velocities and densities
        """
        u_shape = c[('u',0)].shape
        alphaBDF = self.fluidModel.timeIntegration.alpha_bdf
        if  u_shape == self.fluidModel.q[('u',0)].shape:
            vf = self.fluidModel.q[('velocity',0)]
            vs = self.fluidModel.coefficients.q_velocity_solid
            vos = self.fluidModel.coefficients.q_vos
            rho_s = self.fluidModel.coefficients.rho_s
            rho_f = self.fluidModel.coefficients.q_rho
        if  u_shape == self.fluidModel.ebqe[('u',0)].shape:
            vf = self.fluidModel.ebqe[('velocity',0)]
            vs = self.fluidModel.coefficients.ebqe_velocity_solid
            vos = self.fluidModel.coefficients.ebqe_vos
            rho_s = self.fluidModel.coefficients.rho_s
            rho_f = self.fluidModel.coefficients.ebqe_rho
        assert rho_s >= self.rho_s_min, "solid density out of bounds"
        assert (rho_f >= self.rho_f_min).all(), "fluid density out of bounds"
        for i in range(vs.shape[-1]):
            c[('f',0)][...,i] = (1.0-vos)*vf[...,i] + vos*vs[...,i]
        #a is really a scalar diffusion but defining it as diagonal tensor
        #if we push phase momentum interchange (drag) to correction
        #then a may become a full  tensor
        c['a_f'] = 1.0/(self.rho_f_min*alphaBDF)
        c['a_s'] = 1.0/(self.rho_s_min*alphaBDF)
        c[('a',0,0)][...,0] = (1.0-vos)*c['a_f'] + vos*c['a_s']
        for i in range(1,c[('a',0,0)].shape[-1]):
            c[('a',0,0)][...,i] = c[('a',0,0)][...,0]

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

        self.q[('v',0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nDOF_trial_element[0]),
            'd')
        self.q['J'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global,
             self.nSpace_global),
            'd')
        self.q['det(J)'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element),
            'd')
        self.q['inverse(J)'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.nQuadraturePoints_element,
             self.nSpace_global,
             self.nSpace_global),
            'd')
        self.ebq[('v',0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nDOF_trial_element[0]),
            'd')
        self.ebq[('w',0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nDOF_trial_element[0]),
            'd')
        self.ebq['x'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             3),
            'd')
        self.ebq['hat(x)'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             3),
            'd')
        self.ebq['inverse(J)'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global,
             self.nSpace_global),
            'd')
        self.ebq['g'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global-1,
             self.nSpace_global-1),
            'd')
        self.ebq['sqrt(det(g))'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebq['n'] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebq[('dS_u',0)] = numpy.zeros(
            (self.mesh.nElements_global,
             self.mesh.nElementBoundaries_element,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe['dS'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebqe[('dS_u',0)] = self.ebqe['dS']
        self.ebqe['n'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebqe['inverse(J)'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global,
             self.nSpace_global),
            'd')
        self.ebqe['g'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global-1,
             self.nSpace_global-1),
            'd')
        self.ebqe['sqrt(det(g))'] = numpy.zeros(
            (self.mesh.nExteriorElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary),
            'd')
        self.ebq_global['n'] = numpy.zeros(
            (self.mesh.nElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             self.nSpace_global),
            'd')
        self.ebq_global['x'] = numpy.zeros(
            (self.mesh.nElementBoundaries_global,
             self.nElementBoundaryQuadraturePoints_elementBoundary,
             3),
            'd')
        self.ebqe[('advectiveFlux_bc_flag',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('diffusiveFlux_bc_flag',0,0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'i')
        self.ebqe[('advectiveFlux_bc',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux_bc',0,0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('advectiveFlux',0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.ebqe[('diffusiveFlux',0,0)] = numpy.zeros((self.mesh.nExteriorElementBoundaries_global,self.nElementBoundaryQuadraturePoints_elementBoundary),'d')
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
        log("flux bc objects")
        for ci,fbcObject  in self.fluxBoundaryConditionsObjectsDict.iteritems():
            self.ebqe[('advectiveFlux_bc_flag',ci)] = numpy.zeros(self.ebqe[('advectiveFlux_bc',ci)].shape,'i')
            for t,g in fbcObject.advectiveFluxBoundaryConditionsDict.iteritems():
                if self.coefficients.advection.has_key(ci):
                    self.ebqe[('advectiveFlux_bc',ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('advectiveFlux_bc_flag',ci)][t[0],t[1]] = 1
            for ck,diffusiveFluxBoundaryConditionsDict in fbcObject.diffusiveFluxBoundaryConditionsDictDict.iteritems():
                self.ebqe[('diffusiveFlux_bc_flag',ck,ci)] = numpy.zeros(self.ebqe[('diffusiveFlux_bc',ck,ci)].shape,'i')
                for t,g in diffusiveFluxBoundaryConditionsDict.iteritems():
                    self.ebqe[('diffusiveFlux_bc',ck,ci)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
                    self.ebqe[('diffusiveFlux_bc_flag',ck,ci)][t[0],t[1]] = 1
        self.numericalFlux.setDirichletValues(self.ebqe)
        compKernelFlag = 0
        self.presinc = cPresInc.PresInc(
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
        # Load the unknowns into the finite element dof
        self.setUnknowns(u)
        self.numericalFlux.setDirichletValues(self.ebqe)

        # flux boundary conditions
        for t, g in self.fluxBoundaryConditionsObjectsDict[
                0].advectiveFluxBoundaryConditionsDict.iteritems():
            self.ebqe[
                ('advectiveFlux_bc', 0)][
                t[0], t[1]] = g(
                self.ebqe[
                    ('x')][
                    t[0], t[1]], self.timeIntegration.t)
            self.ebqe[('advectiveFlux_bc_flag', 0)][t[0], t[1]] = 1
        for t,g in self.fluxBoundaryConditionsObjectsDict[
                0].diffusiveFluxBoundaryConditionsDictDict[0].iteritems():
            self.ebqe[('diffusiveFlux_bc',0,0)][t[0],t[1]] = g(self.ebqe[('x')][t[0],t[1]],self.timeIntegration.t)
            self.ebqe[('diffusiveFlux_bc_flag',0,0)][t[0],t[1]] = 1

        self.presinc.calculateResidual(  # element
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
            self.numericalFlux.isDOFBoundary[0],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.u[0].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.coefficients.fluidModel.timeIntegration.alpha_bdf,
            self.coefficients.fluidModel.q[('velocity',0)],
            self.coefficients.fluidModel.coefficients.q_velocity_solid,
            self.coefficients.fluidModel.coefficients.q_vos,
            self.coefficients.fluidModel.coefficients.rho_s,
            self.coefficients.fluidModel.coefficients.q_rho,
            self.coefficients.rho_s_min,
            self.coefficients.rho_f_min,
            self.coefficients.fluidModel.ebqe[('velocity',0)],
            self.coefficients.fluidModel.coefficients.ebqe_velocity_solid,
            self.coefficients.fluidModel.coefficients.ebqe_vos,
            self.coefficients.fluidModel.coefficients.ebqe_rho,
            self.q[('u', 0)],
            self.q[('grad(u)', 0)],
            self.ebqe[('u', 0)],
            self.ebqe[('grad(u)', 0)],
            self.numericalFlux.ebqe[('u',0)],
            self.ebqe[('advectiveFlux', 0)],
            self.ebqe[('diffusiveFlux', 0, 0)],
            self.ebqe[('advectiveFlux_bc', 0)],
            self.offset[0], self.stride[0],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray)
        log("Global residual", level=9, data=r)
        #turn this on to view the global conservation residual
        #it should be the same as the linear solver residual tolerance
        #self.coefficients.massConservationError = fabs(
        #    globalSum(r[:self.mesh.nNodes_owned].sum()))
        #log("   Mass Conservation Error", level=3,
        #    data=self.coefficients.massConservationError)
        self.nonlinear_function_evaluations += 1

    def getJacobian(self, jacobian):
        cfemIntegrals.zeroJacobian_CSR(self.nNonzerosInJacobian, jacobian)
        self.presinc.calculateJacobian(  # element
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
            self.numericalFlux.isDOFBoundary[0],
            self.ebqe[('advectiveFlux_bc_flag', 0)],
            self.u[0].femSpace.dofMap.l2g,
            self.u[0].dof,
            self.coefficients.fluidModel.timeIntegration.alpha_bdf,
            self.coefficients.fluidModel.q[('velocity',0)],
            self.coefficients.fluidModel.coefficients.q_velocity_solid,
            self.coefficients.fluidModel.coefficients.q_vos,
            self.coefficients.fluidModel.coefficients.rho_s,
            self.coefficients.fluidModel.coefficients.q_rho,
            self.coefficients.rho_s_min,
            self.coefficients.rho_f_min,
            self.coefficients.fluidModel.ebqe[('velocity',0)],
            self.coefficients.fluidModel.coefficients.ebqe_velocity_solid,
            self.coefficients.fluidModel.coefficients.ebqe_vos,
            self.coefficients.fluidModel.coefficients.ebqe_rho,
            self.csrRowIndeces[(0, 0)], self.csrColumnOffsets[(0, 0)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray,
            self.csrColumnOffsets_eb[(0, 0)])
        log("Jacobian ", level=10, data=jacobian)
        # mwf decide if this is reasonable for solver statistics
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
        log("initalizing ebqe vectors for post-procesing velocity")
        self.u[0].femSpace.elementMaps.getValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                    self.ebqe['x'])
        self.u[0].femSpace.elementMaps.getJacobianValuesGlobalExteriorTrace(self.elementBoundaryQuadraturePoints,
                                                                            self.ebqe['inverse(J)'],
                                                                            self.ebqe['g'],
                                                                            self.ebqe['sqrt(det(g))'],
                                                                            self.ebqe['n'])
        cfemIntegrals.calculateIntegrationWeights(self.ebqe['sqrt(det(g))'],
                                                  self.elementBoundaryQuadratureWeights[('u',0)],
                                                  self.ebqe[('dS_u',0)])
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        log("setting flux boundary conditions")
        self.fluxBoundaryConditionsObjectsDict = dict([(cj,FluxBoundaryConditions(self.mesh,
                                                                                  self.nElementBoundaryQuadraturePoints_elementBoundary,
                                                                                  self.ebqe[('x')],
                                                                                  self.advectiveFluxBoundaryConditionsSetterDict[cj],
                                                                                  self.diffusiveFluxBoundaryConditionsSetterDictDict[cj]))
                                                       for cj in self.advectiveFluxBoundaryConditionsSetterDict.keys()])
        log("initializing coefficients ebqe")
        self.coefficients.initializeGlobalExteriorElementBoundaryQuadrature(self.timeIntegration.t,self.ebqe)

    def estimate_mt(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def updateAfterMeshMotion(self):
        pass
