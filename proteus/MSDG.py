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
        Transport.OneLevelTransport.calculateElementJacobian(self,skipMassTerms):
        self.transfer_lhs = self.elementJacobian.copy()
        self.transfer_source = []
        for ci in self.coefficients.reaction.keys():
            self.transfer_source[ci] = self.elementResidual[ci].copy()
            self.transfer_source[ci][:] = 0.0
            cfemIntegrals.updateReaction_weak(self.q[('r',ci)],
                                              self.q[('w*dV_r',ci)],
                                              self.transfer_source[ci])
        self.transfer_fluxJacobian = self.fluxJacobian.copy()
        self.transfer_fluxJacobian[:] = 0.0
        self.transfer_fluxJacobian_eb = self.fluxJacobian_eb.copy()
        self.transfer_fluxJacobian_eb[:] = 0.0
        self.transfer_fluxJacobian_hj = self.fluxJacobian_hj.copy()
        self.transfer_fluxJacobian_hj[:] = 0.0
        self.numericalFlux.updateInteriorNumericalFluxJacobianTransfer(self.l2g,self.q,self.ebq,self.ebq_global,self.dphi,
                                                                       self.transfer_fluxJacobian,self.transfer_fluxJacobian_eb,self.transfer_fluxJacobian_hj)
        self.transfer_fluxJacobian_exterior = self.fluxJacobian_exterior.copy()
        self.transfer_fluxJacobian_exterior[:] = 0.0
        self.numericalFlux.updateExteriorNumericalFluxJacobianTransfer(self.l2g,self.q,self.ebq,self.ebq_global,self.dphi,
                                                                       self.transfer_fluxJacobian_exterior,self.transfer_fluxJacobian_eb,self.transfer_fluxJacobian_hj)
        #now "assemble" the flux Jacobians into the transfer matrices
        #then invert to get the transfer operator for the element matrices
        #we would then need to assemble the element Jacobians  into a global Jacobian based on the continuous element maps
