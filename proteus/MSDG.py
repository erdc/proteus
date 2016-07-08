import Transport
import pdb
import numpy as np
from Profiling import logEvent


class OneLevelMSDG(Transport.OneLevelTransport):

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
                 stressFluxBoundaryConditionsSetterDict=None,
                 stabilization=None,
                 shockCapturing=None,
                 conservativeFluxDict=None,
                 numericalFluxType=None,
                 TimeIntegrationClass=None,
                 massLumping=False,
                 reactionLumping=False,
                 options=None,
                 name='defaultName',
                 reuse_trial_and_test_quadrature=False,
                 sd=True,
                 movingDomain=False):  # ,
        Transport.OneLevelTransport.__init__(self,
                                             uDict,
                                             phiDict,
                                             testSpaceDict,
                                             matType,
                                             dofBoundaryConditionsDict,
                                             dofBoundaryConditionsSetterDict,
                                             coefficients,
                                             elementQuadrature,
                                             elementBoundaryQuadrature,
                                             fluxBoundaryConditionsDict,
                                             advectiveFluxBoundaryConditionsSetterDict,
                                             diffusiveFluxBoundaryConditionsSetterDictDict,
                                             stressFluxBoundaryConditionsSetterDict,
                                             stabilization,
                                             shockCapturing,
                                             conservativeFluxDict,
                                             numericalFluxType,
                                             TimeIntegrationClass,
                                             massLumping,
                                             reactionLumping,
                                             options,
                                             name,
                                             reuse_trial_and_test_quadrature,
                                             sd,
                                             movingDomain)
        self.cg_spaces = {}
        self.u_cg = {}
        self.nd = coefficients.nd
        from proteus import FemTools
        for i in range(self.nd + 1):
            self.cg_spaces[i] = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(self.mesh,
                                                                                self.nd) # needs to be changed for more general choice of space
            self.u_cg[i] = FemTools.FiniteElementFunction(
                self.cg_spaces[i], name=self.u[i].name)
        # self.cg_transport = Transport.OneLevelTransport(self.u_cg,
        #                                                 self.u_cg,
        #                                                 self.cg_spaces,
        #                                                 matType,
        #                                                 dofBoundaryConditionsDict,
        #                                                 dofBoundaryConditionsSetterDict,
        #                                                 coefficients,
        #                                                 elementQuadrature,
        #                                                 elementBoundaryQuadrature,
        #                                                 fluxBoundaryConditionsDict,
        #                                                 advectiveFluxBoundaryConditionsSetterDict,
        #                                                 diffusiveFluxBoundaryConditionsSetterDictDict,
        #                                                 stressFluxBoundaryConditionsSetterDict,
        #                                                 stabilization,
        #                                                 shockCapturing,
        #                                                 conservativeFluxDict,
        #                                                 numericalFluxType,
        #                                                 TimeIntegrationClass,
        #                                                 massLumping,
        #                                                 reactionLumping,
        #                                                 options,
        #                                                 name,
        #                                                 reuse_trial_and_test_quadrature,
        #                                                 sd,
        #                                                 movingDomain)
        #self.dg_u_all = np.zeros((self.dim,),'d')
        #self.dg_r_all = np.zeros((self.dim,),'d')
        #self.dim = self.cg_transport.dim
        #self.dim_proc = self.cg_transport.dim_proc
        self.T = []
        self.Trhs = []
        for eN in range(self.mesh.nElements_global):
            self.T.append(np.identity(self.nDOF_trial_element[0] * self.nc, 'd'))
            self.Trhs.append(np.zeros(self.nDOF_test_element[0] * self.nc, 'd'))
        self.I = None
        self.RHS = None

    # this function actually does not touched
    def CG_to_DG(self, cg_dof, dg_dof):
        raise SystemExit(0)
        import sys
        print "get to mult"
        sys.exit()
        logEvent("CG -> DG")
        dg_dof[:] = 0.0
        pdb.set_trace()
        for eN in range(self.mesh.nElements_global):
            for i in range(self.nDOF_trial_element[0]):
                for j in range(self.nDOF_trial_element[0]):
                    # 3 indicates the nc
                    dg_dof[3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0] = reduce(lambda x, y: x + y, map(lambda x, y: x * y, self.T[eN][0 * self.nDOF_test_element[0] + i, range(3) * self.nDOF_trial_element[0] + j],
                                                                                                          [cg_dof[3 * self.u_cg[k].femSpace.dofMap.l2g[eN, j] + k] for k in range(0, 3)]))
                    # dg_dof[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 0]  += (self.T[eN][0*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j]*cg_dof[3*self.u_cg[0].femSpace.dofMap.l2g[eN,j]+0] +
                    #                                                        self.T[eN][0*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j]*cg_dof[3*self.u_cg[1].femSpace.dofMap.l2g[eN,j]+1] +
                    # self.T[eN][0*self.nDOF_test_element[0] + i,
                    # 2*self.nDOF_trial_element[0]+j]*cg_dof[3*self.u_cg[2].femSpace.dofMap.l2g[eN,j]+2])
                    dg_dof[3 * self.u[1].femSpace.dofMap.l2g[eN, i] + 1] += (self.T[eN][1 * self.nDOF_test_element[0] + i, 0 * self.nDOF_trial_element[0] + j] * cg_dof[3 * self.u_cg[0].femSpace.dofMap.l2g[eN, j] + 0] +
                                                                             self.T[eN][1 * self.nDOF_test_element[0] + i, 1 * self.nDOF_trial_element[0] + j] * cg_dof[3 * self.u_cg[1].femSpace.dofMap.l2g[eN, j] + 1] +
                                                                             self.T[eN][1 * self.nDOF_test_element[0] + i, 2 * self.nDOF_trial_element[0] + j] * cg_dof[3 * self.u_cg[2].femSpace.dofMap.l2g[eN, j] + 2])
                    dg_dof[3 * self.u[2].femSpace.dofMap.l2g[eN, i] + 2] += (self.T[eN][2 * self.nDOF_test_element[0] + i, 0 * self.nDOF_trial_element[0] + j] * cg_dof[3 * self.u_cg[0].femSpace.dofMap.l2g[eN, j] + 0] +
                                                                             self.T[eN][2 * self.nDOF_test_element[0] + i, 1 * self.nDOF_trial_element[0] + j] * cg_dof[3 * self.u_cg[1].femSpace.dofMap.l2g[eN, j] + 1] +
                                                                             self.T[eN][2 * self.nDOF_test_element[0] + i, 2 * self.nDOF_trial_element[0] + j] * cg_dof[3 * self.u_cg[2].femSpace.dofMap.l2g[eN, j] + 2])

    # this function actually does not touched
    def DG_to_CG(self, dg_dof, cg_dof):
        print "DG -> CG"
        cg_dof[:] = 0.0
        for eN in range(self.mesh.nElements_global):
            for i in range(self.nDOF_trial_element[0]):
                for j in range(self.nDOF_trial_element[0]):
                    cg_dof[3 * self.u_cg[0].femSpace.dofMap.l2g[eN, i] + 0] += (self.T[eN][0 * self.nDOF_trial_element[0] + j, 0 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[0].femSpace.dofMap.l2g[eN, j] + 0] +
                                                                                self.T[eN][1 * self.nDOF_trial_element[0] + j, 0 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[1].femSpace.dofMap.l2g[eN, j] + 1] +
                                                                                self.T[eN][2 * self.nDOF_trial_element[0] + j, 0 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[2].femSpace.dofMap.l2g[eN, j] + 2])
                    cg_dof[3 * self.u_cg[1].femSpace.dofMap.l2g[eN, i] + 1] += (self.T[eN][0 * self.nDOF_trial_element[0] + j, 1 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[0].femSpace.dofMap.l2g[eN, j] + 0] +
                                                                                self.T[eN][1 * self.nDOF_trial_element[0] + j, 1 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[1].femSpace.dofMap.l2g[eN, j] + 1] +
                                                                                self.T[eN][2 * self.nDOF_trial_element[0] + j, 1 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[2].femSpace.dofMap.l2g[eN, j] + 2])
                    cg_dof[3 * self.u_cg[2].femSpace.dofMap.l2g[eN, i] + 2] += (self.T[eN][0 * self.nDOF_trial_element[0] + j, 2 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[0].femSpace.dofMap.l2g[eN, j] + 0] +
                                                                                self.T[eN][1 * self.nDOF_trial_element[0] + j, 2 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[1].femSpace.dofMap.l2g[eN, j] + 1] +
                                                                                self.T[eN][2 * self.nDOF_trial_element[0] + j, 2 * self.nDOF_test_element[0] + i] * dg_dof[3 * self.u[2].femSpace.dofMap.l2g[eN, j] + 2])

    def calculateElementJacobian(self, skipMassTerms=False):

        import numpy as np
        from numpy import linalg
        import cfemIntegrals
        self.til_M = []
        self.bar_M = []
        self.T = []
        self.Trhs = []
        Transport.OneLevelTransport.calculateElementJacobian(
            self, skipMassTerms)
        self.transfer_lhs = {}
        self.transfer_rhs = {}
        # pdb.set_trace()
        print self.elementJacobian.keys()
        for ci in self.elementJacobian.keys():
            for cj in self.elementJacobian[ci].keys():
                self.transfer_lhs[(ci, cj)] = self.elementJacobian[
                    ci][cj].copy()
                self.transfer_rhs[(ci, cj)] = self.elementJacobian[
                    ci][cj].copy()
                self.transfer_rhs[(ci, cj)][:] = 0.0
        # add two more extra terms, not general
        self.transfer_lhs[(1, 2)] = self.elementJacobian[0][0].copy()
        self.transfer_lhs[(1, 2)][:] = 0.0
        self.transfer_lhs[(2, 1)] = self.elementJacobian[0][0].copy()
        self.transfer_lhs[(2, 1)][:] = 0.0
        #
        self.transfer_rhs[(1, 2)] = self.elementJacobian[0][0].copy()
        self.transfer_rhs[(1, 2)][:] = 0.0
        self.transfer_rhs[(2, 1)] = self.elementJacobian[0][0].copy()
        self.transfer_rhs[(2, 1)][:] = 0.0
        #
        ####

        self.transfer_source = {}
        for ci in self.coefficients.reaction.keys():
            self.transfer_source[ci] = self.elementResidual[ci].copy()
            self.transfer_source[ci][:] = 0.0
            cfemIntegrals.updateReaction_weak(self.q[('r', ci)],
                                              self.q[('w*dV_r', ci)],
                                              self.transfer_source[ci])

        self.rhs_source = []

        # add boundary integral terms
        self.alpha_p = 0.9
        self.alpha_beta = 1.0
        for eN in range(self.mesh.nElements_global):
            for ebN in range(self.mesh.nElementBoundaries_element):
                # assumes equal  order across  components
                for i in range(self.nDOF_test_element[0]):
                    # assumes equal  order across components
                    for j in range(self.nDOF_trial_element[0]):
                        for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                            #                             pdb.set_trace()
                            #                             h =1.0
                            #                             lambda_norm = 1.0
                            h = self.q['abs(det(J))'][eN, k] / \
                                self.ebq['sqrt(det(g))'][eN, ebN, k]
                            lambda_norm = 1.0 / \
                                abs(self.ebq[('dr', 1, 1)][eN, ebN, k])

                            #
                            # lhs matrix
                            #
                            # til_D_k
                            self.transfer_lhs[(0, 0)][eN, i, j] += self.alpha_beta * lambda_norm / h * self.ebq[
                                ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # til_C_k
                            self.transfer_lhs[(0, 1)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                                ('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            self.transfer_lhs[(0, 2)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                                ('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # no til_Bk terms from boundary
                            #
                            # til_Ak
                            self.transfer_lhs[(1, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                                'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            self.transfer_lhs[(2, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                                'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]
                            # add more terms
                            self.transfer_lhs[(1, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                                'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            self.transfer_lhs[(2, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                                'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]

                            #
                            # rhs matrix
                            #
                            # bar_D_k
                            self.transfer_rhs[(0, 0)][eN, i, j] += self.alpha_beta * lambda_norm / h * self.ebq[
                                ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # bar_C_k
                            self.transfer_rhs[(0, 1)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                                ('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            self.transfer_rhs[(0, 2)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                                ('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # bar_Bk
                            self.transfer_rhs[(1, 0)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                                ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            self.transfer_rhs[(2, 0)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                                ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]
                            # bar_Ak
                            self.transfer_rhs[(1, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                                'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            self.transfer_rhs[(2, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                                'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]
                            # add more terms
                            self.transfer_rhs[(1, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                                'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            self.transfer_rhs[(2, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                                'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]

            self.til_M.append(
                np.zeros((3 * self.nDOF_test_element[0], 3 * self.nDOF_trial_element[0]), 'd'))
            self.bar_M.append(
                np.zeros((3 * self.nDOF_test_element[0], 3 * self.nDOF_trial_element[0]), 'd'))
            self.T.append(
                np.zeros((3 * self.nDOF_test_element[0], 3 * self.nDOF_trial_element[0]), 'd'))
            self.rhs_source.append(
                np.zeros((3 * self.nDOF_test_element[0], 1), 'd'))
            self.Trhs.append(
                np.zeros((3 * self.nDOF_test_element[0], ), 'd'))

            for i in range(self.nDOF_test_element[0]):
                for j in range(self.nDOF_trial_element[0]):

                    for ic in range(self.nc):
                        for jc in range(self.nc):
                            # til_M
                            self.til_M[eN][ic * self.nDOF_test_element[0] + i,
                                           jc * self.nDOF_trial_element[0] + j] \
                                = self.transfer_lhs[(ic, jc)][eN, i, j]
                            # bar_M
                            self.bar_M[eN][ic * self.nDOF_test_element[0] + i,
                                           jc * self.nDOF_trial_element[0] + j] \
                                = self.transfer_rhs[(ic, jc)][eN, i, j]

                        # add the elementwise source term in transfer operator
                        self.rhs_source[eN][ic * self.nDOF_test_element[0] + i] = self.transfer_source[ic][eN][i]
                    # T
            # pdb.set_trace()
            self.T[eN] = np.dot(linalg.inv(self.til_M[eN]), self.bar_M[eN])
            self.Trhs[eN] = np.dot(linalg.inv(self.til_M[eN]),self.rhs_source[eN])

        # we would then need to assemble the element Jacobians  into a global Jacobian based on the continuous element maps
        # this approach would need formal MSDG space with knowledge of  coarse and fine  spaces
        # let's first check  if M and T are right

    def getRestriction(self):
        class TransferShell:

            def __init__(self, pde):
                self.pde = pde

            def create(self, A):
                pass

            def mult(self, A, x, y):
                raise SystemExit(0)
                import sys
                print "get to mult"
                sys.exit()
                self.pde.DG_to_CG(x, y)
                from petsc4py import PETSc
        return TransferShell(self)

    def getInterpolation(self):
        # pdb.set_trace()
        from petsc4py import PETSc

        class TransferShell:

            def __init__(self, pde):
                self.pde = pde

            def create(self, A):
                pass

            def mult(self, A, x, y):
                self.pde.CG_to_DG(x, y)
        if not self.I:
            self.I = PETSc.Mat().create()
            self.I.setSizes([[3 * self.u[0].femSpace.dim, 3 * self.u[0].femSpace.dim], \
                             [3 * self.u_cg[0].femSpace.dim, 3 * self.u_cg[0].femSpace.dim]])
            #INPY = np.zeros((3*self.u[0].femSpace.dim,3*self.u_cg[0].femSpace.dim),'d')
            self.I.setType(PETSc.Mat.Type.SEQAIJ)
            self.I.setPreallocationNNZ(1)
            self.I.setOption(
                PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        else:
            self.I.zeroEntries()
            self.I.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True)
        # setup righthandside
        # pdb.set_trace()

        if not self.RHS:
            self.RHS = PETSc.Vec().createWithArray(np.zeros((3 * self.u[0].femSpace.dim,), 'd'))
            # self.RHS = PETSc.Vec().create()
            # self.RHS.setSizes(3 * self.u_cg[0].femSpace.dim,)


        logEvent("Assemble global transfer operator in PETSc")
        for eN in range(self.mesh.nElements_global):
            for i in range(self.nDOF_trial_element[0]): #better to change to test
                for j in range(self.nDOF_trial_element[0]):
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0,
                                    3 *
                                    self.u_cg[0].femSpace.dofMap.l2g[
                                        eN, j] + 0,
                                    self.T[eN][0 * self.nDOF_test_element[0] + i,
                                               0 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0,
                                    3 *
                                    self.u_cg[1].femSpace.dofMap.l2g[
                                        eN, j] + 1,
                                    self.T[eN][0 * self.nDOF_test_element[0] + i,
                                               1 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0,
                                    3 *
                                    self.u_cg[2].femSpace.dofMap.l2g[
                                        eN, j] + 2,
                                    self.T[eN][0 * self.nDOF_test_element[0] + i,
                                               2 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1,
                                    3 *
                                    self.u_cg[0].femSpace.dofMap.l2g[
                                        eN, j] + 0,
                                    self.T[eN][1 * self.nDOF_test_element[0] + i,
                                               0 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1,
                                    3 *
                                    self.u_cg[1].femSpace.dofMap.l2g[
                                        eN, j] + 1,
                                    self.T[eN][1 * self.nDOF_test_element[0] + i,
                                               1 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1,
                                    3 *
                                    self.u_cg[2].femSpace.dofMap.l2g[
                                        eN, j] + 2,
                                    self.T[eN][1 * self.nDOF_test_element[0] + i,
                                               2 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2,
                                    3 *
                                    self.u_cg[0].femSpace.dofMap.l2g[
                                        eN, j] + 0,
                                    self.T[eN][2 * self.nDOF_test_element[0] + i,
                                               0 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2,
                                    3 *
                                    self.u_cg[1].femSpace.dofMap.l2g[
                                        eN, j] + 1,
                                    self.T[eN][2 * self.nDOF_test_element[0] + i,
                                               1 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2,
                                    3 *
                                    self.u_cg[2].femSpace.dofMap.l2g[
                                        eN, j] + 2,
                                    self.T[eN][2 * self.nDOF_test_element[0] + i,
                                               2 * self.nDOF_trial_element[0] + j],
                                    PETSc.InsertMode.ADD_VALUES)
                    #
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 0,
                    #            3*self.u_cg[0].femSpace.dofMap.l2g[eN,j]+0] += self.T[eN][0*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 0,
                    #            3*self.u_cg[1].femSpace.dofMap.l2g[eN,j]+1] += self.T[eN][0*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 0,
                    #            3*self.u_cg[2].femSpace.dofMap.l2g[eN,j]+2] += self.T[eN][0*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 1,
                    #            3*self.u_cg[0].femSpace.dofMap.l2g[eN,j]+0] += self.T[eN][1*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 1,
                    #            3*self.u_cg[1].femSpace.dofMap.l2g[eN,j]+1] += self.T[eN][1*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 1,
                    #            3*self.u_cg[2].femSpace.dofMap.l2g[eN,j]+2] += self.T[eN][1*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 2,
                    #            3*self.u_cg[0].femSpace.dofMap.l2g[eN,j]+0] += self.T[eN][2*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 2,
                    #            3*self.u_cg[1].femSpace.dofMap.l2g[eN,j]+1] += self.T[eN][2*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j]
                    # INPY[3*self.u[0].femSpace.dofMap.l2g[eN,i] + 2,
                    # 3*self.u_cg[2].femSpace.dofMap.l2g[eN,j]+2] +=
                    # self.T[eN][2*self.nDOF_test_element[0] + i,
                    # 2*self.nDOF_trial_element[0]+j]
        self.I.assemble()
        logEvent("Assemble global transfer operator in PETSc finished")
        # assemble RHS
        # pdb.set_trace()
        logEvent("Assemble Transfered source term in PETSc")
        for eN in range(self.mesh.nElements_global):
            for i in range(self.nDOF_test_element[0]):
                self.RHS.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0, self.Trhs[eN][0 * self.nDOF_test_element[0] + i],
                PETSc.InsertMode.ADD_VALUES)
                self.RHS.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1, self.Trhs[eN][1 * self.nDOF_test_element[0] + i],
                PETSc.InsertMode.ADD_VALUES)
                self.RHS.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2, self.Trhs[eN][2 * self.nDOF_test_element[0] + i],
                PETSc.InsertMode.ADD_VALUES)
        self.RHS.assemble()
        logEvent("Assemble Transfered source term in PETSc finished")
        # pdb.set_trace()


        return TransferShell(self), self.I, self.RHS # , INPY
    # def getResidual(self,u,r):
    #    #self.CG_to_DG(u, self.dg_u_all)
    #    Transport.OneLevelTransport.getResidual(self, self.dg_u_all, self.dg_r_all)
    #    #import pdb
    #    #pdb.set_trace()
    #    #self.DG_to_CG(self.dg_r_all, r)
    # def getJacobian(self,jacobian,skipMassTerms=False):
    #    Transport.OneLevelTransport.getJacobian(self,jacobian,skipMassTerms)
    #    S,P,I = self.getInterpolation()
    #    import pdb
    #    pdb.set_trace()
