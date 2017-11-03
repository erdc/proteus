from proteus import *
from proteus.default_n import *
from ls_consrv_vortex_2d_p import *
from vortex2D import *


timeIntegrator = ForwardIntegrator
timeIntegration = NoIntegration
stepController = MCorr.Newton_controller#need a tricked up controller that can fix the VOF model's initial conditions

if cDegree_ls==0:
    if useHex:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
        elementQuadrature = CubeGaussQuadrature(nd,vortex_quad_order)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,vortex_quad_order)
    else:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
        base_quad_rule = SimplexGaussQuadrature(nd,vortex_quad_order)
        elementQuadrature = CompositeTriangle(base_quad_rule,hk)
        #elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)
    if parallel or LevelModelType in [MCorr.LevelModel]:#,MCorrElement.LevelModel]:
        numericalFluxType = DoNothing#Diffusion_IIPG_exterior
elif cDegree_ls==-1:
    if pDegree_ls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG


# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None

massLumping = False

shockCapturing = None

multilevelNonlinearSolver  = NLNI

if correctionType == 'dg':
    levelNonlinearSolver = MCorr.ElementNewton
elif correctionType == 'dgp0':
    levelNonlinearSolver = MCorr.ElementConstantNewton
elif correctionType == 'global':
    levelNonlinearSolver = MCorr.GlobalConstantNewton
elif correctionType == 'none':
    levelNonlinearSolver = MCorr.DummyNewton
else:
    levelNonlinearSolver = Newton
    nonlinearSolverNorm = MCorr.conservationNorm
nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolConservation
useEisenstatWalker = True

maxNonlinearIts = 100

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py#PETSc
    levelLinearSolver = KSP_petsc4py#PETSc
    linear_solver_options_prefix = 'mcorr_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU


conservativeFlux = {}
if checkMass:
    auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryMC("vortex2d"+`lRefinement`+"p"+`pDegree_ls`)]
