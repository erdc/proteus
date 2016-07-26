import Transport
import pdb
import numpy as np
from Profiling import logEvent
import LinearAlgebraTools


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
        # add the choice of the cg space
        # import pdb; pdb.set_trace()
        if 'TransferredSpace' in dir(options):
            logEvent('Setting cg space')
            for i in range(self.nd + 1):
                self.cg_spaces[i] = options.TransferredSpace[i](self.mesh,self.nd)
                self.u_cg[i] = FemTools.FiniteElementFunction(
                    self.cg_spaces[i], name=self.u[i].name)

        else:
            for i in range(self.nd + 1):
                self.cg_spaces[i] = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(self.mesh,
                                                                                    self.nd) # needs to be changed for more general choice of space
                self.u_cg[i] = FemTools.FiniteElementFunction(
                    self.cg_spaces[i], name=self.u[i].name)

        # local penalty coeffcients setup
        if 'alpha_p' in dir(options):
            self.alpha_p = options.alpha_p
            logEvent('Setting input local alpha_p = {}'.format(self.alpha_p))
        else:
            self.alpha_p = 0.9
            logEvent('Setting default local alpha_p = {}'.format(self.alpha_p))
        if 'alpha_beta' in dir(options):
            self.alpha_beta = options.alpha_beta
            logEvent('Setting input local alpha_beta = {}'.format(self.alpha_beta))
        else:
            self.alpha_beta = 1.0
            logEvent('Setting default local alpha_beta = {}'.format(self.alpha_beta))
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
        logEvent('Initializing MSDG local transfer operators T from global problem')
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
        for ci in self.elementJacobian.keys(): # equation number
            for cj in self.elementJacobian[ci].keys(): # component number
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

        #specify the size of the mapping matrices
        self.til_M = [np.zeros((self.nc * self.nDOF_test_element[0], self.nc * self.nDOF_trial_element[0]), 'd') for iM in range(self.mesh.nElements_global)]
        self.bar_M = [np.zeros((self.nc * self.nDOF_test_element[0], self.nc * self.cg_spaces[0].max_nDOF_element), 'd') for iM in range(self.mesh.nElements_global)]
        self.T = [np.zeros((self.nc * self.nDOF_test_element[0], self.nc * self.cg_spaces[0].max_nDOF_element), 'd') for iM in range(self.mesh.nElements_global)]
        self.rhs_source = [np.zeros((self.nc * self.nDOF_test_element[0], 1), 'd') for iM in range(self.mesh.nElements_global)]
        self.Trhs = [np.zeros((self.nc * self.nDOF_test_element[0], ), 'd') for iM in range(self.mesh.nElements_global)]

        logEvent('Initializing MSDG local transfer operator T finished')
        logEvent('Setting the extra boundary terms in MSDG local Transfer operator T')
        # add boundary integral terms
        # self.alpha_p = 0.9
        # self.alpha_beta = 1.0

        # debug use
        # logEvent("Initilize debug partiale terms")
        # self.debugD = [np.zeros((6,6), 'd'),np.zeros((6,6), 'd')]
        # self.debugC = [np.zeros((6,12), 'd'),np.zeros((6,12), 'd')]
        # self.debugB = [np.zeros((12,6), 'd'),np.zeros((12,6), 'd')]
        # self.debugA = [np.zeros((12,12), 'd'),np.zeros((12,12), 'd')]
        # self.debugTAss = [np.zeros((18, 18), 'd'), np.zeros((18, 18), 'd')]
        # self.debugTtotal = np.zeros((36,27), 'd')
        # debug ends
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
                            for ic in range(self.nd):
                                self.transfer_lhs[(0, ic + 1)][eN, i, j] -= self.ebq['n'][eN, ebN, k, ic] * self.ebq[
                                    ('v', ic + 1)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]

                            # explicit expression of til_C_k
                            # self.transfer_lhs[(0, 1)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                            #     ('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # self.transfer_lhs[(0, 2)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                            #     ('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]

                            # # no til_Bk terms from boundary
                            #
                            # til_Ak
                            # eps = 1e-5
                            for ia in range(self.nd):
                                for ja in range(self.nd):
                                    # change the penalty term
                                    # if self.ebq['n'][eN, ebN, k, 0] < eps or self.ebq['n'][eN, ebN, k, ia] < eps:
                                    #     self.transfer_lhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm \
                                    #         * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia +1)][eN, ebN, k, i]
                                    # else:
                                    #     self.transfer_lhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, ia] * self.ebq[
                                    #         'n'][eN, ebN, k, ja] * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia +1)][eN, ebN, k, i]





                                    # self.transfer_lhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm \
                                    #     * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia +1)][eN, ebN, k, i]

                                    self.transfer_lhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, ia] * self.ebq[
                                        'n'][eN, ebN, k, ja] * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia +1)][eN, ebN, k, i]
                                    #
                            # explicit til_Ak terms
                            # self.transfer_lhs[(1, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                            #     'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            # self.transfer_lhs[(2, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                            #     'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]

                            # # add more terms

                            # self.transfer_lhs[(1, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                            #     'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            # self.transfer_lhs[(2, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                            #     'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]

                            #
                            # rhs matrix
                            #
                            # bar_D_k
                            self.transfer_rhs[(0, 0)][eN, i, j] += self.alpha_beta * lambda_norm / h * self.ebq[
                                ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # bar_C_k
                            for ic in range(self.nd):
                                    self.transfer_rhs[(0, ic + 1)][eN, i, j] -= self.ebq['n'][eN, ebN, k, ic] * self.ebq[
                                        ('v', ic + 1)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            #explicit expression of bar_C_k
                            # self.transfer_rhs[(0, 1)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                            #     ('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # self.transfer_rhs[(0, 2)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                            #     ('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_H', 0)][eN, ebN, k, i]
                            # bar_Bk
                            for ib in range(self.nd):
                                    self.transfer_rhs[(ib + 1, 0)][eN, i, j] -= self.ebq['n'][eN, ebN, k, ib] * self.ebq[
                                        ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_f', ib + 1)][eN, ebN, k, i]
                            # #explicit expression for bar_BK
                            # self.transfer_rhs[(1, 0)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                            #     ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            # self.transfer_rhs[(2, 0)][eN, i, j] -= self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                            #     ('v', 0)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]

                            # bar_Ak
                            for ia in range(self.nd):
                                for ja in range(self.nd):
                                    # if self.ebq['n'][eN, ebN, k, 0] < eps or self.ebq['n'][eN, ebN, k, ia] < eps:
                                    #     self.transfer_rhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm \
                                    #         * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia + 1)][eN, ebN, k, i]
                                    # else:
                                    #     self.transfer_rhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, ia] * self.ebq[
                                    #         'n'][eN, ebN, k, ja] * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia + 1)][eN, ebN, k, i]





                                    # self.transfer_rhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm \
                                    #     * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia + 1)][eN, ebN, k, i]

                                    self.transfer_rhs[(ia + 1, ja + 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, ia] * self.ebq[
                                        'n'][eN, ebN, k, ja] * self.ebq[('v', ja + 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', ia + 1)][eN, ebN, k, i]
                                    #
                            # # explicit expression for bar_Ak
                            # self.transfer_rhs[(1, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                            #     'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            # self.transfer_rhs[(2, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                            #     'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]
                            # # add more terms
                            # self.transfer_rhs[(1, 2)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 0] * self.ebq[
                            #     'n'][eN, ebN, k, 1] * self.ebq[('v', 2)][eN, ebN, k, j] * self.ebq[('w*dS_f', 1)][eN, ebN, k, i]
                            # self.transfer_rhs[(2, 1)][eN, i, j] += self.alpha_p * h / lambda_norm * self.ebq['n'][eN, ebN, k, 1] * self.ebq[
                            #     'n'][eN, ebN, k, 0] * self.ebq[('v', 1)][eN, ebN, k, j] * self.ebq[('w*dS_f', 2)][eN, ebN, k, i]

            # self.til_M.append(
            #     np.zeros((3 * self.nDOF_test_element[0], 3 * self.nDOF_trial_element[0]), 'd'))
            # self.bar_M.append(
            #     np.zeros((3 * self.nDOF_test_element[0], 3 * self.nDOF_trial_element[0]), 'd'))
            # self.T.append(
            #     np.zeros((3 * self.nDOF_test_element[0], 3 * self.nDOF_trial_element[0]), 'd'))
            # self.rhs_source.append(
            #     np.zeros((3 * self.nDOF_test_element[0], 1), 'd'))
            # self.Trhs.append(
            #     np.zeros((3 * self.nDOF_test_element[0], ), 'd'))

            for i in range(self.nDOF_test_element[0]):
                for j in range(self.nDOF_trial_element[0]):

                    for ic in range(self.nc):
                        for jc in range(self.nc):
                            # til_M
                            self.til_M[eN][self.nc * i + ic,
                                           self.nc * j + jc] \
                                = self.transfer_lhs[(ic, jc)][eN, i, j]
                            # self.til_M[eN][ic * self.nDOF_test_element[0] + i,
                            #                jc * self.nDOF_trial_element[0] + j] \
                            #     = self.transfer_lhs[(ic, jc)][eN, i, j]
                            #
                            # # bar_M
                            # self.bar_M[eN][ic * self.nDOF_test_element[0] + i,
                            #                jc * self.nDOF_trial_element[0] + j] \
                            #     = self.transfer_rhs[(ic, jc)][eN, i, j]

                        # # add the elementwise source term in transfer operator
                        # self.rhs_source[eN][ic * self.nDOF_test_element[0] + i] = self.transfer_source[ic][eN][i]
                    # T
            # pdb.set_trace()
            for i in range(self.nDOF_test_element[0]):
                for j in range(self.cg_spaces[0].max_nDOF_element):
                    for ic in range(self.nc):
                        for jc in range(self.nc):
                            # bar_M
                            self.bar_M[eN][self.nc * i + ic,
                                           self.nc * j + jc] \
                                = self.transfer_rhs[(ic, jc)][eN, i, j]
                            # self.bar_M[eN][ic * self.nDOF_test_element[0] + i,
                            #                jc * self.cg_spaces[0].max_nDOF_element + j] \
                            #     = self.transfer_rhs[(ic, jc)][eN, i, j]

            # debug use
            # for i in range(self.nDOF_test_element[0]):
            #     for j in range(self.cg_spaces[0].max_nDOF_element):
            #         self.debugD[eN][i, j] = self.transfer_rhs[(0, 0)][eN, i, j]
            #         self.debugC[eN][i, j] = self.transfer_rhs[(0, 1)][eN, i, j]
            #         self.debugC[eN][i, 6 + j] = self.transfer_rhs[(0, 2)][eN, i, j]
            #         self.debugB[eN][i, j] = self.transfer_rhs[(1, 0)][eN, i, j]
            #         self.debugB[eN][6 + i, j] = self.transfer_rhs[(2, 0)][eN, i, j]
            #         self.debugA[eN][i, j] = self.transfer_rhs[(1, 1)][eN, i, j]
            #         self.debugA[eN][i, 6 + j] = self.transfer_rhs[(1, 2)][eN, i, j]
            #         self.debugA[eN][6 + i, j] = self.transfer_rhs[(2, 1)][eN, i, j]
            #         self.debugA[eN][6 + i, 6 + j] = self.transfer_rhs[(2, 2)][eN, i, j]
            # self.debugTAss[eN] = np.bmat([[self.debugA[eN], self.debugB[eN]], [self.debugC[eN], self.debugD[eN]]])

            #debug ends

            # add the elementwise source term in transfer operator
            for i in range(self.nDOF_test_element[0]):

                for ic in range(self.nc):
                    self.rhs_source[eN][self.nc * i + ic] = self.transfer_source[ic][eN][i]
                    # self.rhs_source[eN][ic * self.nDOF_test_element[0] + i] = self.transfer_source[ic][eN][i]

            self.T[eN] = np.dot(linalg.inv(self.til_M[eN]), self.bar_M[eN])
            self.Trhs[eN] = np.dot(linalg.inv(self.til_M[eN]),self.rhs_source[eN])
            # import pdb; pdb.set_trace()
        #debug use
        # for eN in range(self.mesh.nElements_global):
        #     for i in range(self.nDOF_trial_element[0]): #better to change to test
        #         for j in range(self.cg_spaces[0].max_nDOF_element):
        #             # D matrix
        #             self.debugTtotal[24 + self.u[0].femSpace.dofMap.l2g[eN,i], 18 + self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugD[eN][i,j]
        #             # C matrix
        #             self.debugTtotal[24 + self.u[0].femSpace.dofMap.l2g[eN,i],  self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugC[eN][i,j]
        #             self.debugTtotal[24 + self.u[0].femSpace.dofMap.l2g[eN,i],  9 + self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugC[eN][i,j+6]
        #
        #             # B matrix
        #
        #             self.debugTtotal[self.u[0].femSpace.dofMap.l2g[eN,i],  18 + self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugB[eN][i,j]
        #             self.debugTtotal[12 + self.u[0].femSpace.dofMap.l2g[eN,i], 18 + self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugB[eN][i + 6,j]
        #
        #             # A matrix
        #             self.debugTtotal[self.u[0].femSpace.dofMap.l2g[eN,i],  self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugA[eN][i,j]
        #             self.debugTtotal[12 +  self.u[0].femSpace.dofMap.l2g[eN,i],  self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugA[eN][i + 6,j]
        #             self.debugTtotal[ self.u[0].femSpace.dofMap.l2g[eN,i],  9 + self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugA[eN][i,j + 6]
        #             self.debugTtotal[12 + self.u[0].femSpace.dofMap.l2g[eN,i],  9 + self.u_cg[0].femSpace.dofMap.l2g[eN, j]] = \
        #             self.debugA[eN][i+6,j+6]
        #

        # debug ends
        logEvent('Setting the extra boundary terms in MSDG local Transfer operator T finished')
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
        # The transfer shell never get used
        class TransferShell:

            def __init__(self, pde):
                self.pde = pde

            def create(self, A):
                pass

            def mult(self, A, x, y):
                self.pde.CG_to_DG(x, y)
        if not self.I:
            self.I = PETSc.Mat().create()
            # import pdb; pdb.set_trace()
            self.I.setSizes([[self.nc * self.u[0].femSpace.dim, self.nc * self.u[0].femSpace.dim], \
                             [self.nc * self.u_cg[0].femSpace.dim, self.nc * self.u_cg[0].femSpace.dim]])
            # self.I.setSizes([[3 * self.u[0].femSpace.dim, 3 * self.u[0].femSpace.dim], \
            #                  [3 * self.u_cg[0].femSpace.dim, 3 * self.u_cg[0].femSpace.dim]])
            #
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
            self.RHS = PETSc.Vec().createWithArray(np.zeros((self.nc * self.u[0].femSpace.dim,), 'd'))
            # self.RHS = PETSc.Vec().createWithArray(np.zeros((3 * self.u[0].femSpace.dim,), 'd'))
            # self.RHS = PETSc.Vec().create()
            # self.RHS.setSizes(3 * self.u_cg[0].femSpace.dim,)

        logEvent("Assemble global transfer operator in PETSc")
        # import pdb; pdb.set_trace()
        self.debugT = np.zeros((self.nc * self.u[0].femSpace.dim, self.nc * self.u_cg[0].femSpace.dim), 'd')
        for eN in range(self.mesh.nElements_global):
            for i in range(self.nDOF_trial_element[0]): #better to change to test
                for j in range(self.cg_spaces[0].max_nDOF_element):
                    #implicit expression of assembly
                    for ic in range(self.nc):
                        for jc in range(self.nc):
                            self.I.setValue(self.nc * self.u[0].femSpace.dofMap.l2g[eN, i] + ic,
                                            self.nc * self.u_cg[0].femSpace.dofMap.l2g[eN, j] + jc,
                                            self.T[eN][self.nc * i + ic,
                                                       self.nc * j + jc],
                                            PETSc.InsertMode.ADD_VALUES)
                            # import pdb; pdb.set_trace()
                            self.debugT[self.nc * self.u[0].femSpace.dofMap.l2g[eN, i] + ic,
                                        self.nc * self.u_cg[0].femSpace.dofMap.l2g[eN, j] + jc]  = \
                                        self.T[eN][self.nc * i + ic,
                                                   self.nc * j + jc]
                    # explicit expression of assembly
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0,
                    #                 3 *
                    #                 self.u_cg[0].femSpace.dofMap.l2g[
                    #                     eN, j] + 0,
                    #                 self.T[eN][0 * self.nDOF_test_element[0] + i,
                    #                            0 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0,
                    #                 3 *
                    #                 self.u_cg[1].femSpace.dofMap.l2g[
                    #                     eN, j] + 1,
                    #                 self.T[eN][0 * self.nDOF_test_element[0] + i,
                    #                            1 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0,
                    #                 3 *
                    #                 self.u_cg[2].femSpace.dofMap.l2g[
                    #                     eN, j] + 2,
                    #                 self.T[eN][0 * self.nDOF_test_element[0] + i,
                    #                            2 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1,
                    #                 3 *
                    #                 self.u_cg[0].femSpace.dofMap.l2g[
                    #                     eN, j] + 0,
                    #                 self.T[eN][1 * self.nDOF_test_element[0] + i,
                    #                            0 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1,
                    #                 3 *
                    #                 self.u_cg[1].femSpace.dofMap.l2g[
                    #                     eN, j] + 1,
                    #                 self.T[eN][1 * self.nDOF_test_element[0] + i,
                    #                            1 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1,
                    #                 3 *
                    #                 self.u_cg[2].femSpace.dofMap.l2g[
                    #                     eN, j] + 2,
                    #                 self.T[eN][1 * self.nDOF_test_element[0] + i,
                    #                            2 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2,
                    #                 3 *
                    #                 self.u_cg[0].femSpace.dofMap.l2g[
                    #                     eN, j] + 0,
                    #                 self.T[eN][2 * self.nDOF_test_element[0] + i,
                    #                            0 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2,
                    #                 3 *
                    #                 self.u_cg[1].femSpace.dofMap.l2g[
                    #                     eN, j] + 1,
                    #                 self.T[eN][2 * self.nDOF_test_element[0] + i,
                    #                            1 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)
                    # self.I.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2,
                    #                 3 *
                    #                 self.u_cg[2].femSpace.dofMap.l2g[
                    #                     eN, j] + 2,
                    #                 self.T[eN][2 * self.nDOF_test_element[0] + i,
                    #                            2 * self.nDOF_trial_element[0] + j],
                    #                 PETSc.InsertMode.ADD_VALUES)



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
        # self.debugT_PETSc = LinearAlgebraTools.petsc4py_sparse_2_dense(self.I, output = 'transferMatrix')
        logEvent("Assemble global transfer operator in PETSc finished")
        # assemble RHS
        # pdb.set_trace()
        logEvent("Assemble Transfered source term in PETSc")
        for eN in range(self.mesh.nElements_global):
            for i in range(self.nDOF_test_element[0]):

                for ic in range(self.nc):
                    self.RHS.setValue(self.nc * self.u[0].femSpace.dofMap.l2g[eN, i] + ic,
                                      self.Trhs[eN][self.nc * i + ic],
                                      PETSc.InsertMode.ADD_VALUES)
                # # explicit assemly RHS
                # self.RHS.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 0, self.Trhs[eN][0 * self.nDOF_test_element[0] + i],
                # PETSc.InsertMode.ADD_VALUES)
                # self.RHS.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 1, self.Trhs[eN][1 * self.nDOF_test_element[0] + i],
                # PETSc.InsertMode.ADD_VALUES)
                # self.RHS.setValue(3 * self.u[0].femSpace.dofMap.l2g[eN, i] + 2, self.Trhs[eN][2 * self.nDOF_test_element[0] + i],
                # PETSc.InsertMode.ADD_VALUES)
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
