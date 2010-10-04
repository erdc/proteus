"""
A class hierarchy for shock capturing diffusion methods
"""
import numpy
import cshockCapturing
class ShockCapturing_base:
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True):
        self.nc = coefficients.nc
        self.nd = nd
        self.components=range(self.nc)
        self.shockCapturingFactor=shockCapturingFactor
        self.lag=lag
        #mwf had to put in da_sge, df_sge
        self.coefficients = coefficients
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff=[]
        self.numDiff_last=[]
        for ci in range(self.nc):
            if self.lag:
                self.numDiff_last.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(numpy.zeros(cq[('u',ci)].shape,'d'))
            else:
                self.numDiff_last.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(cq[('numDiff',ci,ci)])
        #mwf what if just want shock capturing and no stabilization?
        #need to put in lag steps for da_sge and df_sge below as well if not
        #done by Subgrid error
        for ci in range(self.nc):
            for cj in range(self.nc):
                if cq.has_key(('df',ci,cj)) and not cq.has_key(('df_sge',ci,cj)):
                    cq[('df_sge',ci,cj)]=cq[('df',ci,cj)]
        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                if not cq.has_key(('grad(phi)_sge',ck)):
                    cq[('grad(phi)_sge',ck)]=cq[('grad(phi)',ck)]
                for cj in cjDict.keys():
                    if not cq.has_key(('dphi_sge',ck,cj)):
                        cq[('dphi_sge',ck,cj)]=cq[('dphi',ck,cj)]
                    if not cq.has_key(('da_sge',ci,ck,cj)):
                        cq[('da_sge',ci,ck,cj)]=cq[('da',ci,ck,cj)]
        
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            self.numDiff[ci].flat[:] = self.shockCapturingFactor
    def updateShockCapturingHistory(self):
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]

ConstantDiffusion_SC = ShockCapturing_base

class ResGrad_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusionResGrad(self.shockCapturingFactor,
                                                               self.mesh.elementDiametersArray,
                                                               q[('pdeResidual',ci)],
                                                               q[('grad(u)',ci)],
                                                               self.numDiff[ci])
class ResGradFFDarcy_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
    def calculateNumericalDiffusion(self,q):
        cshockCapturing.calculateNumericalDiffusionResGrad(self.shockCapturingFactor,
                                                           self.mesh.elementDiametersArray,
                                                           q[('pdeResidual',0)],
                                                           q[('grad(u)',0)],
                                                           self.numDiff[0])
        
class ResGradQuad_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,gradLag=True):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.debug=False
	self.gradLag = gradLag
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            if self.debug:
                if numpy.isnan(q[('pdeResidual',ci)]).any(): # 
                    import pdb
                    print "NaN's in res"
                    pdb.set_trace()
            cshockCapturing.calculateNumericalDiffusionResGradQuad(self.shockCapturingFactor,
                                                                   self.mesh.elementDiametersArray,
                                                                   q[('pdeResidual',ci)],
                                                                   q[('grad(u)',ci)],
                                                                   self.numDiff[ci])
            if self.debug:
                if numpy.isnan(self.numDiff[ci]).any():
                    import pdb
                    print "NaN's in numDiff"
                    pdb.set_trace()
    #todo ido, fix this so it doesn't break models that don't have a 'gradNorm'
    # def initializeElementQuadrature(self,mesh,t,cq):
    #     ShockCapturing_base.initializeElementQuadrature(self,mesh,t,cq)
    #     self.gradNorm=[]
    #     self.gradNorm_last=[]
    #     for ci in range(self.nc):
    #         if self.gradLag:
    #             self.gradNorm_last.append(cq[('gradNorm',ci,ci)])
    #             self.gradNorm.append(numpy.zeros(cq[('u',ci)].shape,'d'))
    #         else:
    #             self.gradNorm_last.append(cq[('gradNorm',ci,ci)])
    #             self.gradNorm.append(cq[('gradNorm',ci,ci)])

    # def updateShockCapturingHistory(self):
    #     ShockCapturing_base.updateShockCapturingHistory(self)
    #     if self.gradLag:
    #         for ci in range(self.nc):
    #             self.gradNorm_last[ci][:] = self.gradNorm[ci]

class Eikonal_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.debug=False
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusionEikonal(self.shockCapturingFactor,
                                                               self.mesh.elementDiametersArray,
                                                               q[('pdeResidual',ci)],
                                                               self.numDiff[ci])

class ScalarAdvection_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.5,lag=True):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusion_A_1(self.shockCapturingFactor,
                                                            self.mesh.elementDiametersArray,
                                                            q[('pdeResidual',ci)],
                                                            q[('mt',ci)],
                                                            q[('df',ci,ci)],
                                                            self.numDiff[ci])

class HamiltonJacobi_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFlag='1',shockCapturingFactor=0.5,lag=True):
        self.shockCapturingFlag=shockCapturingFlag=shockCapturingFlag
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusionHJ(self.shockCapturingFlag,
                                                          self.shockCapturingFactor,
                                                          self.mesh.elementDiametersArray,
                                                          q[('pdeResidual',ci)],
                                                          q[('mt',ci)],
                                                          q[('H',ci)],
                                                          self.numDiff[ci])
class HamiltonJacobiJaffre_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFlag='1',shockCapturingFactor=0.5,lag=True,betaPower=0.1):
        self.shockCapturingFlag=shockCapturingFlag=shockCapturingFlag
        self.beta=betaPower
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusionJaffre(self.shockCapturingFactor,
                                                              self.beta,
                                                              self.mesh.elementDiametersArray,
                                                              q[('pdeResidual',ci)],
                                                              q[('dH',ci,ci)],
                                                              self.numDiff[ci])
class JaffreGradU_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFlag='1',shockCapturingFactor=0.5,lag=True,betaPower=0.1):
        self.shockCapturingFlag=shockCapturingFlag=shockCapturingFlag
        self.beta=betaPower
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusionJaffre(self.shockCapturingFactor,
                                                              self.beta,
                                                              self.mesh.elementDiametersArray,
                                                              q[('pdeResidual',ci)],
                                                              q[('grad(u)',ci)],
                                                              self.numDiff[ci])

class ResGradJuanes_SC(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,uSC=1.0,lag=True):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.uSC = uSC
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusionResGradJuanes(self.shockCapturingFactor,
                                                                     self.uSC,
                                                                     self.mesh.elementDiametersArray,
                                                                     q[('pdeResidual',ci)],
                                                                     q[('grad(u)',ci)],
                                                                     self.numDiff[ci])

#examples of shock capturing with lagging allowed after certain number of steps
class ResGradDelayLag_SC(ResGrad_SC):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,nStepsToDelay=None):
        ResGrad_SC.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff=[]
        self.numDiff_last=[]
        self.cq_numDiff=[]
        for ci in range(self.nc):
            if self.lag:
                self.numDiff_last.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(numpy.zeros(cq[('u',ci)].shape,'d'))
            elif self.lag == False and self.nStepsToDelay != None:
                self.cq_numDiff.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(cq[('numDiff',ci,ci)])
            else:
                self.numDiff.append(cq[('numDiff',ci,ci)])
        #mwf what if just want shock capturing and no stabilization?
        #need to put in lag steps for da_sge and df_sge below as well if not
        #done by Subgrid error
        for ci in range(self.nc):
            for cj in range(self.nc):
                if cq.has_key(('df',ci,cj)) and not cq.has_key(('df_sge',ci,cj)):
                    cq[('df_sge',ci,cj)]=cq[('df',ci,cj)]
        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                if not cq.has_key(('grad(phi)_sge',ck)):
                    cq[('grad(phi)_sge',ck)]=cq[('grad(phi)',ck)]
                for cj in cjDict.keys():
                    if not cq.has_key(('dphi_sge',ck,cj)):
                        cq[('dphi_sge',ck,cj)]=cq[('dphi',ck,cj)]
                    if not cq.has_key(('da_sge',ci,ck,cj)):
                        cq[('da_sge',ci,ck,cj)]=cq[('da',ci,ck,cj)]
    #
    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay != None and self.nSteps > self.nStepsToDelay:
            self.lag = True
            self.numDiff = []
            self.numDiff_last=[]
            for ci in range(self.nc):
                self.numDiff_last.append(self.cq_numDiff[ci])
                self.numDiff.append(numpy.zeros(self.cq_numDiff[ci].shape,'d'))

class ResGradQuadDelayLag_SC(ResGradQuad_SC):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,nStepsToDelay=None):
        ResGradQuad_SC.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff=[]
        self.numDiff_last=[]
        self.cq_numDiff=[]
        for ci in range(self.nc):
            if self.lag:
                self.numDiff_last.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(numpy.zeros(cq[('u',ci)].shape,'d'))
            elif self.lag == False and self.nStepsToDelay != None:
                self.cq_numDiff.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(cq[('numDiff',ci,ci)])
            else:
                self.numDiff.append(cq[('numDiff',ci,ci)])
    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay != None and self.nSteps > self.nStepsToDelay:
            self.lag = True
            self.numDiff = []
            self.numDiff_last=[]
            for ci in range(self.nc):
                self.numDiff_last.append(self.cq_numDiff[ci])
                self.numDiff.append(numpy.zeros(self.cq_numDiff[ci].shape,'d'))

class NavierStokes_SC(ResGradQuad_SC):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,nStepsToDelay=None):
        ResGradQuad_SC.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
    def calculateNumericalDiffusion(self,q):
        for ci in range(1,self.nc):
            if numpy.isnan(q[('pdeResidual',ci)]).any():
                import pdb
                print "NaN's in res"
                pdb.set_trace()
            cshockCapturing.calculateNumericalDiffusionResGradQuad(self.shockCapturingFactor,
                                                                   self.mesh.elementDiametersArray,
                                                                   q[('pdeResidual',ci)],
                                                                   q[('grad(u)',ci)],
                                                                   self.numDiff[ci])
            if numpy.isnan(self.numDiff[ci]).any():
                import pdb
                print "NaN's in numDiff"
                pdb.set_trace()
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff={}
        self.numDiff_last={}
        self.cq_numDiff={}
        for ci in range(1,self.nc):
            if self.lag:
                self.numDiff_last[ci]=cq[('numDiff',ci,ci)]
                self.numDiff[ci]=numpy.zeros(cq[('u',ci)].shape,'d')
            elif self.lag == False and self.nStepsToDelay != None:
                self.cq_numDiff[ci]=cq[('numDiff',ci,ci)]
                self.numDiff[ci]=cq[('numDiff',ci,ci)]
            else:
                self.numDiff[ci]=cq[('numDiff',ci,ci)]
    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(1,self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay != None and self.nSteps > self.nStepsToDelay:
            self.lag = True
            self.numDiff = []
            self.numDiff_last=[]
            for ci in range(1,self.nc):
                self.numDiff_last[ci]=self.cq_numDiff[ci]
                self.numDiff[ci]=numpy.zeros(self.cq_numDiff[ci].shape,'d')

class NavierStokes_SC_opt(ResGradQuad_SC):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,nStepsToDelay=None):
        ResGradQuad_SC.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
    def calculateNumericalDiffusion(self,q):
        for ci in range(1,self.nc):
            if numpy.isnan(q[('pdeResidual',ci)]).any():
                import pdb
                print "NaN's in res"
                pdb.set_trace()
            cshockCapturing.calculateNumericalDiffusionResGradQuad(self.shockCapturingFactor,
                                                                   self.mesh.elementDiametersArray,
                                                                   q[('pdeResidual',ci)],
                                                                   q[('grad(u)',ci)],
                                                                   self.numDiff[ci])
            if numpy.isnan(self.numDiff[ci]).any():
                import pdb
                print "NaN's in numDiff"
                pdb.set_trace()
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff={}
        self.numDiff_last={}
        for ci in range(1,self.nc):
            if self.lag:
                self.numDiff_last[ci]=numpy.zeros(cq[('u',ci)].shape,'d')
                self.numDiff[ci]=cq[('numDiff',ci,ci)]
            else:
                self.numDiff_last[ci]=cq[('numDiff',ci,ci)]
                self.numDiff[ci]=cq[('numDiff',ci,ci)]
    def updateShockCapturingHistory(self):
        for ci in range(1,self.nc):
            if self.lag:
                self.numDiff_last[ci][:]=self.numDiff[ci]

class ResGradFFDarcyDelayLag_SC(ResGradFFDarcy_SC):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,nStepsToDelay=None):
        ResGradFFDarcy_SC.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.nStepsToDelay = nStepsToDelay
        self.nSteps=0
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff=[]
        self.numDiff_last=[]
        self.cq_numDiff=[]
        for ci in range(self.nc):
            if self.lag:
                self.numDiff_last.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(numpy.zeros(cq[('u',ci)].shape,'d'))
            elif self.lag == False and self.nStepsToDelay != None:
                self.cq_numDiff.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(cq[('numDiff',ci,ci)])
            else:
                self.numDiff.append(cq[('numDiff',ci,ci)])
    def updateShockCapturingHistory(self):
        self.nSteps += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        if self.lag == False and self.nStepsToDelay != None and self.nSteps > self.nStepsToDelay:
            self.lag = True
            self.numDiff = []
            self.numDiff_last=[]
            for ci in range(self.nc):
                self.numDiff_last.append(self.cq_numDiff[ci])
                self.numDiff.append(numpy.zeros(self.cq_numDiff[ci].shape,'d'))
                

