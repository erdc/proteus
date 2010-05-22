from pyadh import *
from pyadh.default_n import *
from ls_vortex_3d_p import *
from vortex import *
nd = 3

if timeIntegration_ls == "BE":
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    #timeIntegration = VBDF
    #stepController = Min_dt_cfl_controller
    #timeOrder =2
elif timeIntegration_ls == "FLCBDF":
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller
elif timeIntegration_ls == "RK":
    if cDegree_ls == -1:
        timeIntegration = LinearSSPRKPIintegration
    else:
        timeIntegration = LinearSSPRKintegration        
    stepController=Min_dt_RKcontroller
    timeOrder = pDegree_ls+1
    nStagesTime = timeOrder
else:
    raise RuntimeError

if pDegree_ls == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elif pDegree_ls == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
else:
    print "pDegree_ls = %s not recognized " % pDegree_ls
#femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)


subgridError = None
if useHJ:
    if LevelModelType == NCLS.OneLevelNCLS:
        subgridError = HamiltonJacobi_ASGS_opt(coefficients,nd,lag=True)
    else:
        subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)
else:
    subgridError = Advection_ASGS(coefficients,nd,lag=True)

massLumping = False

numericalFluxType = None

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=shockCapturingFactor_ls,lag=True)
if LevelModelType == NCLS.OneLevelNCLS or parallel:
        numericalFluxType = DoNothing


multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 0.01/(nn -1.0 )

maxNonlinearIts = 50

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU


linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {}

#checkMass = True

if not applyCorrection:
   auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryLS("vortex3dnc"+`lRefinement`)]
