from proteus import *
from proteus.default_n import *
from MCorr_p import *
from cons_ls import *

multilevelNonlinearSolver  = NLNI
if ct.STABILIZATION_TYPE_vof==0: #SUPG methods 
    levelNonlinearSolver = Newton
else: #Edge based stabilization
    levelNonlinearSolver = NewtonWithL2ProjectionForMassCorrection

nonlinearSolverNorm = MCorr.conservationNorm
nonlinearSmoother = NLGaussSeidel
fullNewtonFlag = True

timeIntegrator = ForwardIntegrator
timeIntegration = NoIntegration
stepController = MCorr.Newton_controller#need a tricked up controller that can fix the VOF model's initial conditions

if useHex:
    if pDegree_ncls==1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_ncls==2:
        femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    elementQuadrature = CubeGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,quad_order)
else:
    if pDegree_ncls==1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ncls==2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}        
    elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)
if parallel or LevelModelType in [MCorr.LevelModel]:#,MCorrElement.LevelModel]:
    numericalFluxType = DoNothing#Diffusion_IIPG_exterior

subgridError = None
massLumping = False
shockCapturing = None

tolFac = 0.0
nl_atol_res = atolConservation
useEisenstatWalker = True

maxNonlinearIts = 100

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'mcorr_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU
    
conservativeFlux = {}
if checkMass:
    auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryMC("vortex2d"+`lRefinement`+"p"+`pDegree_ls`)]
