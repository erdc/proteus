from pyadh import *
from pyadh.default_n import *
from redist_square_1d_p import *
from square import *


timeIntegration = NoIntegration
timeIntegrator = SteadyStateIntegrator
stepController = Newton_controller

if LevelModelType == RDLS.OneLevelRDLS or True:
    timeIntegration = BackwardEuler
    stepController  = Osher_PsiTC_controller
    runCFL=1.0
    rtol_res[0] = 0.0
    atol_res[0] = 0.01/float(nn-1)#1.0e-6
    
if cDegree_ls==0:
    if pDegree_ls==1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ls==2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    #original 
    subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=True)
    if LevelModelType == RDLS.OneLevelRDLS or True:
        subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='2',lag=True)#mwf debug HamiltonJacobi_ASGS_opt(coefficients,nd,stabFlag='2',lag=True)
        shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.0,lag=True)
        numericalFluxType = DoNothing
elif cDegree_ls==-1:
    timeIntegration = BackwardEuler_cfl
    stepController = Osher_PsiTC_controller
    runCFL=0.1
    rtol_res[0] = 0.0
    atol_res[0] = 0.01/float(nn-1)#1.0e-6
    if pDegree_ls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    numericalFluxType=HamiltonJacobi_DiagonalChengShu
    #numericalFluxtype=HamiltonJacobi_DiagonalLesaintRaviart
    #shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=shockCapturingFactor_rd,lag=False)
    
elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolRedistance

maxNonlinearIts = 50 #1 for PTC

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {}
