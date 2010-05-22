from pyadh import *
from pyadh.default_n import *
from transport_biodegradation01_2d_p import *

timeIntegration = FLCBDF
stepController = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3
rtol_u[1] = 1.0e-3
atol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
atol_u[2] = 1.0e-3

DT = None
runCFL = 0.9
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

#mwf everything gets SimplexGaussQuadrature
#cek makeing simple 1c problem
elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)
#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)


nLevels = 1
massLumping=False
#subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton
fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU
#multilevelLinearSolver = NI
#numericalFlux = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

levelLinearSolver = LU
#levelLinearSolver = MGM
#levelLinearSolver = StarILU
#levelLinearSolver = GaussSeidel
#levelLinearSolver = Jacobi

smoother = StarILU
smoother = GaussSeidel
smoother = Jacobi

linTolFac = 1.0e-10

conservativeFlux = None
