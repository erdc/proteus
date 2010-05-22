from pyadh import *
from pyadh.default_n import *
from LN_Example1_adr_p import *
#BackwardEuler
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

if useELLAM:
    shockCapturing = None
    subgridError   = None
else:
    subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)
    shockCapturing = ResGrad_SC(coefficients,nd,lag=True)
    
    
multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
maxNonlinearIts = 1

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

atol = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

