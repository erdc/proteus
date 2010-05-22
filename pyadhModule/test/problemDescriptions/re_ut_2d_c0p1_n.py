from pyadh import *
from pyadh.default_n import *
from re_ut_2d_p import *

nnx=21
nny=21
nLevels = 1

DT = 0.001
nDTout = int(T/DT)
timeIntegrator  = ForwardIntegrator
timeIntegration = BackwardEuler
timeIntegration = FLCBDF
stepController = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-2
atol_u[0] = 1.0e-2

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#appears to work
elementQuadrature = {'default':SimplexLobattoQuadrature(nd,3),
                    ('u',0):SimplexGaussQuadrature(nd,3)}
elementBoundaryQuadrature ={'default':SimplexGaussQuadrature(nd-1,2),
                           ('u',0):SimplexGaussQuadrature(nd-1,2)}

#mwf looks like works too
#elementQuadrature = {'default':SimplexGaussQuadrature(nd,3),
#                     ('u',0):SimplexGaussQuadrature(nd,3)}
#elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3),
#                             ('u',0):SimplexGaussQuadrature(nd-1,3)}

class AdvectionDiffusionReaction_ASGSv2(SGE_base):
    def __init__(self,coefficients,nd,stabFlag='1',lag=True,turnOnLag=False,nNoLag=3):
        SGE_base.__init__(self,coefficients,nd,lag)
        self.stabilizationFlag = stabFlag
        self.turnOnLag = turnOnLag
        self.nNoLag    = nNoLag
        self.nUpdates  = 0
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
            
subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = AdvectionDiffusionReaction_ASGSv2(coefficients,nd,stabFlag='2',lag=False,turnOnLag=True)

class ResGrad_SCv2(SC_base):
    def __init__(self,coefficients,nd,shockCapturingFactor=0.25,lag=True,turnOnLag=False,nNoLag=3):
        SC_base.__init__(self,coefficients,nd,shockCapturingFactor,lag)
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
    def calculatenumpyalDiffusion(self,q):
        for ci in range(self.nc):
            cshockCapturing.calculatenumpyalDiffusionResGrad(self.shockCapturingFactor,
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
#shockCapturing = ResGrad_SCv2(coefficients,nd,shockCapturingFactor=0.9,lag=False,turnOnLag=True)


massLumping = False

numericalFluxType = None

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
#multilevelNonlinearSolver = Newton

#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

#nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 1.0e-4

nl_atol_res = 1.0e-4

maxNonlinearIts = 10#0#1

maxLineSearches = 100
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

