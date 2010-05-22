from pyadh import *
from pyadh.default_n import *
from beach_erosion_board_waves import *
from vof_beach_erosion_board_waves_2d_p import *


timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)


subgridError = None
subgridError = Advection_ASGS(coefficients,nd,lag=False)

massLumping = False

#numericalFluxType = None
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

shockCapturing = None

#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=True)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atolVolumeOfFluid


maxNonlinearIts = 50

matrix = SparseMatrix

if usePETSc:
    multilevelLinearSolver= PETSc

    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
