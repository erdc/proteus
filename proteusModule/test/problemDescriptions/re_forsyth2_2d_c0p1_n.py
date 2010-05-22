from pyadh import *
from pyadh.default_n import *
from re_forsyth2_2d_p import *

nn=3
nLevels = 1
#level 1 : 97, 257, 161, 2.82843 
#level 2 : 354, 997, 644, 1.41421
#level 3 :1351,3926 ,2576, 0.707107 
#level 4 :5277,  15580, 10304, 0.353553 
#level 5 :20857,  62072, 41216, 0.176777 

hLevels = {0:2.82843,1:1.41421,2:0.707107,3:0.353553,4:0.176777}
triangleOptions="pAq30Dena%f" % (0.5*hLevels[2]**2,)

#timeIntegration = BackwardEuler
#DT = 0.01*hLevels[nLevels-1]/abs(rechargeRate)
#DT = 0.01
#DT = 0.001
nDTout = 100#int(T/DT)
#runCFL=0.1

DT = None#0.001
#timeIntegrator  = ForwardIntegrator
#timeIntegration = BackwardEuler
timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
atol_u[0] = 1.0e-5
rtol_u[0] = 1.0e-5

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#appears to work
# elementQuadrature = {'default':SimplexLobattoQuadrature(nd,3),
#                     ('u',0):SimplexGaussQuadrature(nd,3)}
# elementBoundaryQuadrature ={'default':SimplexGaussQuadrature(nd-1,2),
#                            ('u',0):SimplexGaussQuadrature(nd-1,2)}

#mwf looks like works too
# elementQuadrature = {'default':SimplexGaussQuadrature(nd,3),
#                      ('u',0):SimplexGaussQuadrature(nd,3)}
# elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3),
#                              ('u',0):SimplexGaussQuadrature(nd-1,3)}

elementQuadrature=SimplexLobattoQuadrature(nd,3)
elementBoundaryQuadrature=SimplexLobattoQuadrature(nd-1,3)

class AdvectionDiffusionReaction_ASGSv2(SGE_base):
    def __init__(self,coefficients,nd,stabFlag='1',lag=True,turnOnLag=False,nNoLag=3):
        SGE_base.__init__(self,coefficients,nd,lag)
        self.stabilizationFlag = stabFlag
        self.turnOnLag = turnOnLag
        self.nNoLag    = nNoLag
        self.nUpdates  = 0
        self.coefficients=coefficients
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.tau=[]
        self.tau_last=[]
        for ci in range(self.nc):
            if self.lag or self.turnOnLag:
                self.tau_last.append(numpy.zeros(cq[('u',ci)].shape,'d'))
                self.tau.append(numpy.zeros(cq[('u',ci)].shape,'d'))
            else:
                self.tau.append(numpy.zeros(cq[('u',ci)].shape,'d'))
            if cq.has_key(('df',ci,ci)):
                self.df_last = copy.deepcopy(cq[('df',ci,ci)])
                cq[('df_sge',ci,ci)] = self.df_last
        self.cq=cq
        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                cq[('grad(phi)_sge',ck)]=copy.deepcopy(cq[('grad(phi)',ck)])
                for cj in cjDict.keys():
                    cq[('dphi_sge',ck,cj)]=copy.deepcopy(cq[('dphi',ck,cj)])
                    cq[('da_sge',ci,ck,cj)]=copy.deepcopy(cq[('da',ci,ck,cj)])
    def calculateSubgridError(self,q):
        for ci in range(self.nc):
            csubgridError.calculateSubgridError_ADR_tau(self.stabilizationFlag,
                                                        self.mesh.elementDiametersArray,
                                                        q[('dmt',ci,ci)],
                                                        q[('df',ci,ci)],
                                                        q[('a',ci,ci)],
                                                        q[('da',ci,ci,ci)],
                                                        q[('grad(phi)',ci)],
                                                        q[('dphi',ci,ci)],
                                                        q[('dr',ci,ci)],
                                                        q[('pe',ci)],
                                                        q[('cfl',ci)],
                                                        self.tau[ci])
            if self.lag:
                tau=self.tau_last[ci]
            else:
                q[('df_sge',ci,ci)][:] = q[('df',ci,ci)]
                q[('da_sge',ci,ci,ci)][:]=q[('da',ci,ci,ci)]
                q[('grad(phi)_sge',ci)][:]=q[('grad(phi)',ci)]
                q[('dphi_sge',ci,ci)][:]=q[('dphi',ci,ci)]
                tau=self.tau[ci]
            csubgridError.calculateSubgridError_tauRes(tau,
                                                       q[('pdeResidual',ci)],
                                                       q[('dpdeResidual',ci,ci)],
                                                       q[('subgridError',ci)],
                                                       q[('dsubgridError',ci,ci)])
#             print "tau",tau
#             print "pdeResidual",q[('pdeResidual',ci)]
#             print "dpdeResidual",q[('dpdeResidual',ci,ci)]
#             print "subgrid error",q[('subgridError',ci)]
#             print "dsubgrid error",q[('dsubgridError',ci,ci)]
    def updateSubgridErrorHistory(self):
        self.nUpdates += 1
        if self.lag:
            for ci in range(self.nc):
                self.tau_last[ci][:] = self.tau[ci]
        elif self.turnOnLag == True:
            self.lag = self.nUpdates > self.nNoLag
            for ci in range(self.nc):
                self.tau_last[ci][:] = self.tau[ci]
                self.df_last[:] = self.cq[('df',ci,ci)]
        for ci,ckDict in self.coefficients.diffusion.iteritems():
            for ck,cjDict in ckDict.iteritems():
                self.cq[('grad(phi)_sge',ck)][:]=self.cq[('grad(phi)',ck)]
                for cj in cjDict.keys():
                    self.cq[('dphi_sge',ck,cj)][:]=self.cq[('dphi',ck,cj)]
                    self.cq[('da_sge',ci,ck,cj)][:]=self.cq[('da',ci,ck,cj)]
            
subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = AdvectionDiffusionReaction_ASGSv2(coefficients,nd,stabFlag='2',lag=False,turnOnLag=True)


class ResGrad_SCv2(ShockCapturing_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,turnOnLag=False,nNoLag=3):
        ShockCapturing_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
        self.turnOnLag = turnOnLag
        self.nNoLag   = nNoLag
        self.nUpdates = 0
    def initializeElementQuadrature(self,mesh,t,cq):
        self.mesh=mesh
        self.numDiff=[]
        self.numDiff_last=[]
        for ci in range(self.nc):
            if self.lag or self.turnOnLag:
                self.numDiff_last.append(cq[('numDiff',ci,ci)])
                self.numDiff.append(numpy.zeros(cq[('u',ci)].shape,'d'))
            else:
                self.numDiff.append(cq[('numDiff',ci,ci)])
        self.cq=cq
    def calculateNumericalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculateNumericalDiffusionResGrad(self.shockCapturingFactor,
                                                               self.mesh.elementDiametersArray,
                                                               q[('pdeResidual',ci)],
                                                               q[('grad(u)',ci)],
                                                               self.numDiff[ci])
    def updateShockCapturingHistory(self):
        self.nUpdates += 1
        if self.lag:
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
        elif self.turnOnLag == True:
            self.lag = self.nUpdates > self.nNoLag
            for ci in range(self.nc):
                self.numDiff_last[ci][:] = self.numDiff[ci]
shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)

#shockCapturing = ResGrad_SCv2(coefficients,nd,shockCapturingFactor=0.5,lag=False,turnOnLag=True)
shockCapturing = ResGradDelayLag_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False,nStepsToDelay=10)

massLumping = False

numericalFluxType = None

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
#multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
multilevelNonlinearSolver = Newton

#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

#nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

maxNonlinearIts = 10#0#1

maxLineSearches = 10#0
matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = MGM

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

#conservativeFlux = {0:'sun-gs-rt0'}#{0:'pwl'}#{0:'pwc'}#
#if using sun-gs-rt0 need weak dirichlet conditions
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

