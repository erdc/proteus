
class LAstab1C:
    def __init__(self,mesh,nc,nd,stencil):
        self.mesh = mesh
        self.nc = nc
        self.nd = nd
        self.tau=None
        self.tau_last = None
        self.lag_tau=True
    def calculateSubgridError(self,q):
        ci=0
        if self.tau == None:
            self.tau = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
            self.tau_last = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
            if self.lag_tau:
                vfemIntegrals.calculateSubgridError_A_tau(self.mesh.elementDiametersArray,
                                                          q[('df',ci,ci)],
                                                          q[('dmt',ci,ci)],
                                                          self.tau_last)

        #go ahead and compute cfl
        vfemIntegrals.calculateSubgridError_A_tau(self.mesh.elementDiametersArray,
                                                  q[('df',ci,ci)],
                                                  q[('dmt',ci,ci)],
                                                  self.tau)
        if (self.lag_tau):
            vfemIntegrals.calculateSubgridError_tau_res(self.tau_last,
                                                        q[('pdeResidual',ci)],
                                                        q[('dpdeResidual',ci,ci)],
                                                        q[('subgridError',ci)],
                                                        q[('dsubgridError',ci,ci)])
        else:
            vfemIntegrals.calculateSubgridError_tau_res(self.tau,
                                                        q[('pdeResidual',ci)],
                                                        q[('dpdeResidual',ci,ci)],
                                                        q[('subgridError',ci)],
                                                        q[('dsubgridError',ci,ci)])
        vfemIntegrals.calculateCFLADR(self.mesh.elementDiametersArray,
                                      q[('dm',ci,ci)],
                                      q[('df',ci,ci)],
                                      q[('cfl',ci)]);
    def updateSubgridErrorHistory(self):
        if self.tau != None:
            self.tau_last[:] = self.tau

class LAstab:
    def __init__(self,mesh,nc,nd,stencil,ci=0):
        self.mesh = mesh
        self.nc = nc
        self.nd = nd
        self.tau={}
        self.tau_last = {}
        self.lag_tau=True
        for ci in range(self.nc):
            self.tau[ci]=None; self.tau_last[ci]=None
    def calculateSubgridError(self,q):
        for ci in range(self.nc):
            if self.tau[ci] == None:
                self.tau[ci] = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
                self.tau_last[ci] = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
                if self.lag_tau:
                    vfemIntegrals.calculateSubgridError_A_tau(self.mesh.elementDiametersArray,
                                                              q[('df',ci,ci)],
                                                              q[('dmt',ci,ci)],
                                                              self.tau_last[ci])

            #go ahead and compute cfl
            vfemIntegrals.calculateSubgridError_A_tau(self.mesh.elementDiametersArray,
                                                      q[('df',ci,ci)],
                                                      q[('dmt',ci,ci)],
                                                      self.tau[ci])
            if (self.lag_tau):
                vfemIntegrals.calculateSubgridError_tau_res(self.tau_last[ci],
                                                            q[('pdeResidual',ci)],
                                                            q[('dpdeResidual',ci,ci)],
                                                            q[('subgridError',ci)],
                                                            q[('dsubgridError',ci,ci)])
            else:
                vfemIntegrals.calculateSubgridError_tau_res(self.tau[ci],
                                                            q[('pdeResidual',ci)],
                                                            q[('dpdeResidual',ci,ci)],
                                                            q[('subgridError',ci)],
                                                            q[('dsubgridError',ci,ci)])
            vfemIntegrals.calculateCFLADR(self.mesh.elementDiametersArray,
                                          q[('dm',ci,ci)],
                                          q[('df',ci,ci)],
                                          q[('cfl',ci)])
        #enc ci
    #def
    def updateSubgridErrorHistory(self):
        for ci in range(self.nc):
            if self.tau[ci] != None:
                self.tau_last[ci][:] = self.tau[ci]


class LAshock:
    def __init__(self,mesh,nc,nd,stencil):
        self.mesh = mesh
        self.nc = nc
        self.nd = nd
        self.numDiff={}
        self.numDiff_last = {}
        self.shockCapturingFactor={}
        for ci in range(self.nc):
            self.numDiff[ci]=None; self.numDiff_last[ci]=None;
            self.shockCapturingFactor[ci]=0.5
        #end ci
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            if self.numDiff[ci] == None:
                self.numDiff[ci] = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
                self.numDiff_last[ci] = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
                vfemIntegrals.calculateNumericalDiffusion_A_1(self.shockCapturingFactor[ci],
                                                              self.mesh.elementDiametersArray,
                                                              q[('pdeResidual',ci)],
                                                              q[('mt',ci)],
                                                              q[('df',ci,ci)],
                                                              self.numDiff_last[ci])
            #end if
            vfemIntegrals.calculateNumericalDiffusion_A_1(self.shockCapturingFactor[ci],
                                                              self.mesh.elementDiametersArray,
                                                              q[('pdeResidual',ci)],
                                                              q[('mt',ci)],
                                                              q[('df',ci,ci)],
                                                              self.numDiff[ci])
            q[('numDiff',ci,ci)].flat[:]=self.numDiff_last[ci].flat
        #end for
    #end def
    def updateShockCapturingHistory(self):
        for ci in range(self.nc):
            if self.numDiff[ci] != None:
                self.numDiff_last[ci][:] = self.numDiff[ci]
            #if
        #for
    #def

class LAshock1C:
    def __init__(self,mesh,nc,nd,stencil,ci=0):
        self.mesh = mesh
        self.nc = nc
        self.nd = nd
        self.numDiff=None
        self.numDiff_last = None
        self.shockCapturingFactor=0.5
        self.ci=ci
    def calculateNumericalDiffusion(self,q):
        ci=self.ci
        if self.numDiff == None:
            self.numDiff = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
            self.numDiff_last = Numeric.zeros(q[('dmt',ci,ci)].shape,Numeric.Float)
        vfemIntegrals.calculateNumericalDiffusion_A_1(self.shockCapturingFactor,
                                                      self.mesh.elementDiametersArray,
                                                      q[('pdeResidual',ci)],
                                                      q[('mt',ci)],
                                                      q[('df',ci,ci)],
                                                      self.numDiff)
        q[('numDiff',ci,ci)].flat[:]=self.numDiff_last.flat      
    def updateShockCapturingHistory(self):
        if self.numDiff != None:
            self.numDiff_last[:] = self.numDiff


