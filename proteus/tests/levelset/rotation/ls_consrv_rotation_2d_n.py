from __future__ import absolute_import
from proteus import *
from proteus.default_n import *
try:
    from .ls_consrv_rotation_2d_p import *
    from .rotation2D import *
except:
    from ls_consrv_rotation_2d_p import *
    from rotation2D import *


timeIntegrator = ForwardIntegrator
timeIntegration = NoIntegration
stepController = MCorr.Newton_controller#need a tricked up controller that can fix the VOF model's initial conditions

if cDegree_ls==0:
    if useHex:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
        elementQuadrature = CubeGaussQuadrature(nd,rotation_quad_order)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,rotation_quad_order)
    else:
        if pDegree_ls==1:
            femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
        elif pDegree_ls==2:
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}        
        elementQuadrature = SimplexGaussQuadrature(nd,rotation_quad_order)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,rotation_quad_order)
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
    auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryMC("rotation2d"+repr(lRefinement)+"p"+repr(pDegree_ls))]
