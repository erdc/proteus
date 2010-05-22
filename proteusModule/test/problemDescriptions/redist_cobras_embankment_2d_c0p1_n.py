from pyadh import *
from pyadh.default_n import *
from redist_cobras_embankment_2d_p import *
from embankmentRunup import *

if rdtimeIntegration == 'newton':    
    timeIntegration = NoIntegration
    stepController = Newton_controller
elif rdtimeIntegration == 'tte':
    timeIntegration = BackwardEuler_cfl
    timeIntegration = PsiTCtte
elif rdtimeIntegration == 'osher-fmm':
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_FMM_controller
    runCFL=1.0
else:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_controller
    runCFL=1.0
#     timeIntegration = PsiTCtte
#     stepController = PsiTCtte_controller
#     rtol_res[0] = 0.0
#     atol_res[0] = 0.1*L[0]/(nn-1.0)#10% of he
#runCFL=1.0
#DT=None

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
if spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

subgridError = None
if rdtimeIntegration == 'newton':
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=False)
else:
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=True)
    
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=True)

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=False)
if rdtimeIntegration == 'newton':    
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=False)
else:
    shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=rd_shockCapturingFactor,lag=True)
    
massLumping = False

numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind

multilevelNonlinearSolver  = MultilevelEikonalSolver
levelNonlinearSolver = UnstructuredFMMandFSWsolvers.FMMEikonalSolver
multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton
if rdtimeIntegration != 'newton':    
    maxLineSearches = 0
nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

#this needs to be set appropriately for pseudo-transient
tolFac = 0.0

nl_atol_res = 0.01*L[0]/nn

atol_res[0] = 1.0e-6 #for pseudo transient
rtol_res[0] = 0.0

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

if usePETSc:
    numericalFluxType = NF_base# does nothing

    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
