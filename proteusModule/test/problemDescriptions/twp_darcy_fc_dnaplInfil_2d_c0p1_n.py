from pyadh import *
from pyadh.default_n import *
from twp_darcy_fc_dnaplInfil_2d_p import *

timeIntegrator = ForwardIntegrator
#mwf works ok
#timeIntegration= BackwardEuler
#stepController = Min_dt_controller
#runCFL=0.1
#if useForsyth:
#    DT = 1.0e0#1.0e2#1.0e-2
#else:
#    DT = 10.0
#nDTout = 50#int(T/DT)
#runCFL=None


timeIntegration = FLCBDF_TwophaseDarcy_fc
stepController = FLCBDF_controller
rtol_u[0] = 1.0e-4
atol_u[0] = 1.0e-4
rtol_u[1] = 1.0e-4
atol_u[1] = 1.0e-4
DT = None
nDTout = 50#int(T/DT)


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

elementQuadrature = SimplexLobattoQuadrature(nd,1)
elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nn=41
nLevels = 1
triangleOptions = "q30Dena0.005A"
nLevels = 1

#need to handle nondiagonal parts?
class TwophaseFC_AdvectionDiffusionReaction_ASGS(SGE_base):
    def __init__(self,coefficients,nd,stabFlag='1',lag=True):
        SGE_base.__init__(self,coefficients,nd,lag)
        self.stabilizationFlag = stabFlag
    def calculateSubgridError(self,q):
        ci = 0 #wetting
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
        ci = 1 #wetting
        csubgridError.calculateSubgridError_ADR_tau(self.stabilizationFlag,
                                                    self.mesh.elementDiametersArray,
                                                    q[('dmt',ci,0)],#need wrt saturation?
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


subgridError = None
massLumping=False
shockCapturing = None

#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)


multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 10#025
maxLineSearches = 10#0

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = PETSc#LU
#for petsc do things lie
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol  1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
#can also set -pc_asm_overlap 2 with default asm type (restrict)

levelLinearSolver = PETSc#LU
#pick number of layers to use in overlap 
#nLayersOfOverlapForParallel = 1
#parallelPartitioningType = MeshParallelPartitioningTypes.node
parallelPartitioningType = MeshParallelPartitioningTypes.element
numericalFluxType = DarcyFC_IIPG_exterior

linTolFac = 0.0001

conservativeFlux = {0:'pwl',1:'pwl'}
