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
import cMy_correction


class NumericalFlux(proteus.NumericalFlux.ConstantAdvection_Diffusion_SIPG_exterior):
    def __init__(self,
                 vt,
                 getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        proteus.NumericalFlux.ConstantAdvection_Diffusion_SIPG_exterior.__init__(self, vt, getPointwiseBoundaryConditions,
                                                                                 getAdvectiveFluxBoundaryConditions,
                                                                                 getDiffusiveFluxBoundaryConditions)


class Coefficients(TC_base):
    """
    return divergence free U solving
    .. math::
       u = U+\nabla p
    with pure Nuemann boundary condition
    """

    def __init__(self,
                 nd=2,
                 modelIndex=None,
                 pModelIndex=None,
                 vModelIndex=None):

        assert(nd in [2, 3])
        self.nd = nd
        if self.nd == 2:
            sdInfo = {(0, 0): (np.array([0, 1, 2], dtype='i'),
                               np.array([0, 1], dtype='i'))}
        else:
            sdInfo = {(0, 0): (np.array([0, 1, 2, 3], dtype='i'),
                               np.array([0, 1, 2], dtype='i'))}
        
        mass = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction = {}
        hamiltonian = {}
        if nd == 2:
            variableNames = ['correct_u', 'correct_v']
            advection = {0: {0: 'constant',
                             1: 'constant'},
                         1: {0: 'constant',
                             1: 'constant'}}
            TC_base.__init__(self,
                             2,
                             mass,
                             advection,
                             diffusion,
                             potential,
                             reaction,
                             hamiltonian,
                             variableNames,
                             sparseDiffusionTensors=sdInfo,
                             useSparseDiffusion=True)
            self.vectorComponents = [0, 1]
            self.vectorName = "correct_velocity"

        self.modelIndex = modelIndex
        self.pModelIndex = pModelIndex
        self.vModelIndex = vModelIndex

    def attachModels(self, modelList):
        """
        Attach the model for velocity and density to PresureIncrement model
        """
        self.model = modelList[self.modelIndex]
        self.pModel = modelList[self.pModelIndex]
        self.vModel = modelList[self.vModelIndex]

    def initializeMesh(self, mesh):
        """
        Give the TC object access to the mesh for any mesh-dependent information.
        """
        pass

    def initializeElementQuadrature(self, t, cq):
        """
        Give the TC object access to the element quadrature storage
        """
        pass

    def initializeElementBoundaryQuadrature(self, t, cebq, cebq_global):
        """
        Give the TC object access to the element boundary quadrature storage
        """
        pass

    def initializeGlobalExteriorElementBoundaryQuadrature(self, t, cebqe):
        """
        Give the TC object access to the exterior element boundary quadrature storage
        """
        pass

    def initializeGeneralizedInterpolationPointQuadrature(self, t, cip):
        """
        Give the TC object access to the generalized interpolation point storage. These points are used  to project nonlinear potentials (phi).
        """
        pass

    def preStep(self, t, firstStep=False):
        """
        Move the current values to values_last to keep cached set of values for bdf1 algorithm
        """
        copyInstructions = {}
        return copyInstructions

    def postStep(self, t, firstStep=False):
        """
        Update the fluid velocities
        """
        copyInstructions = {}
        return copyInstructions

    def evaluate(self, t, c):
        self.evaluatePresureIncrement(t, c)

    def evaluatePresureIncrement(self, t, c):
        """
        Evaluate the coefficients after getting the velocities and densities
        """
        import pdb
        pdb.set_trace()
#         u_shape = c[('u', 0)].shape
#         alphaBDF = self.fluidModel.timeIntegration.alpha_bdf
#         if u_shape == self.fluidModel.q[('u', 0)].shape:
#             vf = self.fluidModel.q[('velocity', 0)]
#             vs = self.fluidModel.coefficients.q_velocity_solid
#             vos = self.fluidModel.coefficients.q_vos
#             rho_s = self.fluidModel.coefficients.rho_s
#             rho_f = self.fluidModel.coefficients.q_rho
#         if u_shape == self.fluidModel.ebqe[('u', 0)].shape:
#             vf = self.fluidModel.ebqe[('velocity', 0)]
#             vs = self.fluidModel.coefficients.ebqe_velocity_solid
#             vos = self.fluidModel.coefficients.ebqe_vos
#             rho_s = self.fluidModel.coefficients.rho_s
#             rho_f = self.fluidModel.coefficients.ebqe_rho
# 
#         assert rho_s >= self.rho_s_min, "solid density out of bounds"
#         assert (rho_f >= self.rho_f_min).all(), "fluid density out of bounds"
#         for i in range(vs.shape[-1]):
#             c[('f', 0)][..., i] = (1.0 - vos) * vf[..., i] + vos * vs[..., i]
#         # a is really a scalar diffusion but defining it as diagonal tensor
#         # if we push phase momentum interchange (drag) to correction
#         # then a may become a full  tensor
#         c['a_f'] = 1.0 / (self.rho_f_min * alphaBDF)
#         c['a_s'] = 1.0 / (self.rho_s_min * alphaBDF)
#         c[('a', 0, 0)][..., 0] = (1.0 - vos) * c['a_f'] + vos * c['a_s']
#         for i in range(1, c[('a', 0, 0)].shape[-1]):
#             c[('a', 0, 0)][..., i] = c[('a', 0, 0)][..., 0]


class LevelModel(proteus.Transport.OneLevelTransport):
    nCalls = 0
    
    def __init__(self,*args,**kargs):
        proteus.Transport.OneLevelTransport.__init__(self,*args,**kargs)
        
        compKernelFlag = 0
        self.presinc = cMy_correction.My_correction(
            self.nSpace_global,
            self.nQuadraturePoints_element,
            self.u[0].femSpace.elementMaps.localFunctionSpace.dim,
            self.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
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
            self.coefficients.pModel.u[0].femSpace.psi,
            self.coefficients.pModel.u[0].femSpace.grad_psi,
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
            self.u[0].dof,
            self.u[1].dof,
            self.coefficients.vModel.u[0].dof,
            self.coefficients.vModel.u[1].dof,
            self.coefficients.pModel.u[0].femSpace.dofMap.l2g,
            self.coefficients.pModel.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
            self.coefficients.pModel.u[0].dof,
            self.offset[0], 
            self.stride[0],
            self.offset[1], 
            self.stride[1],
            r,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray)


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
            self.coefficients.pModel.u[0].femSpace.psi,
            self.coefficients.pModel.u[0].femSpace.grad_psi,
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
            self.u[0].dof,
            self.u[1].dof,
            self.coefficients.vModel.u[0].dof,
            self.coefficients.vModel.u[1].dof,
            self.coefficients.pModel.u[0].femSpace.dofMap.l2g,
            self.coefficients.pModel.u[0].femSpace.referenceFiniteElement.localFunctionSpace.dim,
            self.coefficients.pModel.u[0].dof,
            self.csrRowIndeces[(0, 0)], 
            self.csrColumnOffsets[(0, 0)],
            self.csrRowIndeces[(0, 1)], 
            self.csrColumnOffsets[(0, 1)],
            self.csrRowIndeces[(1, 0)], 
            self.csrColumnOffsets[(1, 0)],
            self.csrRowIndeces[(1, 1)], 
            self.csrColumnOffsets[(1, 1)],
            jacobian,
            self.mesh.nExteriorElementBoundaries_global,
            self.mesh.exteriorElementBoundariesArray,
            self.mesh.elementBoundaryElementsArray,
            self.mesh.elementBoundaryLocalElementBoundariesArray)

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
                                                  self.elementBoundaryQuadratureWeights[('u', 0)],
                                                  self.ebqe[('dS_u', 0)])
        self.u[0].femSpace.elementMaps.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.elementMaps.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisValuesTraceRef(
            self.elementBoundaryQuadraturePoints)
        self.u[0].femSpace.getBasisGradientValuesTraceRef(
            self.elementBoundaryQuadraturePoints)

    def estimate_mt(self):
        pass

    def calculateAuxiliaryQuantitiesAfterStep(self):
        pass

    def calculateSolutionAtQuadrature(self):
        pass

    def updateAfterMeshMotion(self):
        pass
