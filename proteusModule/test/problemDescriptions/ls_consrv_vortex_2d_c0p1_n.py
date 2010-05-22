from pyadh import *
from pyadh.default_n import *
from ls_consrv_vortex_2d_p import *
from vortex import *


timeIntegrator = ForwardIntegrator
timeIntegration = NoIntegration
stepController = Newton_controller

if cDegree_ls==0:
    if pDegree_ls==1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ls==2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    if parallel or LevelModelType ==  MCorr.OneLevelMCorr:
        numericalFluxType = DoNothing#Diffusion_IIPG_exterior
elif cDegree_ls==-1:
    if pDegree_ls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG

elementQuadrature = SimplexGaussQuadrature(nd,vortex_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vortex_quad_order)

# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None

massLumping = False

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolConservation

maxNonlinearIts = 100

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU

    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 1.0e-6

conservativeFlux = {}
auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryMC("vortex2d"+`lRefinement`+"p"+`pDegree_ls`)]
