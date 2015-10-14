import  Transport

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
                 sd = True,
                 movingDomain=False):#,
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
    def calculateElementJacobian(self,skipMassTerms=False):
        import numpy as np
        from numpy import  linalg
        import  cfemIntegrals
        self.til_M = []
        self.bar_M = []
        self.T = []
        Transport.OneLevelTransport.calculateElementJacobian(self,skipMassTerms)
        self.transfer_lhs={}
        self.transfer_rhs={}
        for ci in self.elementJacobian.keys():
            for cj in self.elementJacobian[ci].keys():
                self.transfer_lhs[(ci,cj)] = self.elementJacobian[ci][cj].copy()
                self.transfer_rhs[(ci,cj)] = self.elementJacobian[ci][cj].copy()
                self.transfer_rhs[(ci,cj)][:] = 0.0
        self.transfer_source = {}
        for ci in self.coefficients.reaction.keys():
            self.transfer_source[ci] = self.elementResidual[ci].copy()
            self.transfer_source[ci][:] = 0.0
            cfemIntegrals.updateReaction_weak(self.q[('r',ci)],
                                              self.q[('w*dV_r',ci)],
                                              self.transfer_source[ci])
        #add boundary integral terms
        self.alpha_p = 0.5
        self.alpha_beta = 0.5
        for eN in range(self.mesh.nElements_global):
            for  ebN in range(self.mesh.nElementBoundaries_element):
                for i in range(self.nDOF_test_element[0]):#assumes equal  order across  components
                    for j in range(self.nDOF_trial_element[0]):#assumes equal  order across components
                        for k in range(self.nElementBoundaryQuadraturePoints_elementBoundary):
                            #
                            #lhs matrix
                            #
                            #til_D_k
                            self.transfer_lhs[(0,0)][eN,i,j] += self.alpha_beta*self.ebq[('v',0)][eN,ebN,k,j]*self.ebq[('w*dS_H',0)][eN,ebN,k,i]
                            #til_C_k
                            self.transfer_lhs[(0,1)][eN,i,j] -= self.ebq['n'][eN,ebN,k,0]*self.ebq[('v',1)][eN,ebN,k,j]*self.ebq[('w*dS_H',0)][eN,ebN,k,i]
                            self.transfer_lhs[(0,2)][eN,i,j] -= self.ebq['n'][eN,ebN,k,1]*self.ebq[('v',2)][eN,ebN,k,j]*self.ebq[('w*dS_H',0)][eN,ebN,k,i]
                            #no til_Bk terms from boundary
                            #
                            #til_Ak
                            self.transfer_lhs[(1,1)][eN,i,j] += self.alpha_p*self.ebq['n'][eN,ebN,k,0]*self.ebq[('v',1)][eN,ebN,k,j]*self.ebq[('w*dS_f',1)][eN,ebN,k,i]
                            self.transfer_lhs[(2,2)][eN,i,j] += self.alpha_p*self.ebq['n'][eN,ebN,k,1]*self.ebq[('v',2)][eN,ebN,k,j]*self.ebq[('w*dS_f',2)][eN,ebN,k,i]
                            #
                            #rhs matrix
                            #
                            #bar_D_k
                            self.transfer_rhs[(0,0)][eN,i,j] += self.alpha_beta*self.ebq[('v',0)][eN,ebN,k,j]*self.ebq[('w*dS_H',0)][eN,ebN,k,i]
                            #bar_C_k
                            self.transfer_rhs[(0,1)][eN,i,j] -= self.ebq['n'][eN,ebN,k,0]*self.ebq[('v',1)][eN,ebN,k,j]*self.ebq[('w*dS_H',0)][eN,ebN,k,i]
                            self.transfer_rhs[(0,2)][eN,i,j] -= self.ebq['n'][eN,ebN,k,1]*self.ebq[('v',2)][eN,ebN,k,j]*self.ebq[('w*dS_H',0)][eN,ebN,k,i]
                            #bar_Bk
                            self.transfer_rhs[(1,0)] -= self.ebq[('v',0)][eN,ebN,k,j]*self.ebq[('w*dS_f',1)][eN,ebN,k,i]
                            self.transfer_rhs[(2,0)] -= self.ebq[('v',0)][eN,ebN,k,j]*self.ebq[('w*dS_f',2)][eN,ebN,k,i]
                            #bar_Ak
                            self.transfer_rhs[(1,1)][eN,i,j] -= self.alpha_p*self.ebq['n'][eN,ebN,k,0]*self.ebq[('v',1)][eN,ebN,k,j]*self.ebq[('w*dS_f',1)][eN,ebN,k,i]
                            self.transfer_rhs[(2,2)][eN,i,j] -= self.alpha_p*self.ebq['n'][eN,ebN,k,1]*self.ebq[('v',2)][eN,ebN,k,j]*self.ebq[('w*dS_f',2)][eN,ebN,k,i]
            self.til_M.append(np.zeros((3*self.nDOF_test_element[0],3*self.nDOF_trial_element[0]),'d'))
            self.bar_M.append(np.zeros((3*self.nDOF_test_element[0],3*self.nDOF_trial_element[0]),'d'))
            self.T.append(np.zeros((3*self.nDOF_test_element[0],3*self.nDOF_trial_element[0]),'d'))
            for  i in range(self.nDOF_test_element[0]):
                for  j in range(self.nDOF_trial_element[0]):
                    #til_M
                    self.til_M[eN][0*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(0,0)][eN,i,j]
                    self.til_M[eN][0*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(0,1)][eN,i,j]
                    self.til_M[eN][0*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(0,2)][eN,i,j]
                    self.til_M[eN][1*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(1,0)][eN,i,j]
                    self.til_M[eN][1*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(1,1)][eN,i,j]
                    #self.til_M[eN][1*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(1,2)][eN,i,j]
                    self.til_M[eN][2*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(2,0)][eN,i,j]
                    #self.til_M[eN][2*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(2,1)][eN,i,j]
                    self.til_M[eN][2*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j] = self.transfer_lhs[(2,2)][eN,i,j]
                    #bar_M
                    self.bar_M[eN][0*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(0,0)][eN,i,j]
                    self.bar_M[eN][0*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(0,1)][eN,i,j]
                    self.bar_M[eN][0*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(0,2)][eN,i,j]
                    self.bar_M[eN][1*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(1,0)][eN,i,j]
                    self.bar_M[eN][1*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(1,1)][eN,i,j]
                    #self.bar_M[eN][1*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(1,2)][eN,i,j]
                    self.bar_M[eN][2*self.nDOF_test_element[0] + i, 0*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(2,0)][eN,i,j]
                    #self.bar_M[eN][2*self.nDOF_test_element[0] + i, 1*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(2,1)][eN,i,j]
                    self.bar_M[eN][2*self.nDOF_test_element[0] + i, 2*self.nDOF_trial_element[0]+j] = self.transfer_rhs[(2,2)][eN,i,j]
                    #T
            self.T[eN] = linalg.inv(self.til_M[eN])*self.bar_M[eN]
        #we would then need to assemble the element Jacobians  into a global Jacobian based on the continuous element maps
        #this approach would need formal MSDG space with knowledge of  coarse and fine  spaces
        #let's first check  if M and T are right
