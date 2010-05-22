from pyadh import *
from pyadh.default_n import *
from stokes_pome_2d_homog_p import *

timeIntegration = NoIntegration
## timeIntegration = FLCBDF
## stepController  = FLCBDF_controller
## systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
## rtol_u[1] = 1.0e-2
## rtol_u[2] = 1.0e-2
## atol_u[1] = 1.0e-3
## atol_u[2] = 1.0e-3
tnList=[0.0,T]

triangleOptions="YpAq30Dena%f" % (0.25*(dx)**2,)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

#elementQuadrature = SimplexGaussQuadrature(nd,3)

#elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,3)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,3)

#femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#             1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#             2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,5)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

nn=3
nLevels = 1

#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)
subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False

shockCapturing = None



multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-10

matrix = SparseMatrix
#matrix = numpy.array

#PETSc eg #-pc_type lu -pc_factor_mat_solver_package
multilevelLinearSolver = PETSc#LU
multilevelLinearSolver = LU

levelLinearSolver = PETSc#LU
levelLinearSolver = LU
#amount of overlap for parallel 
nLayersOfOverlapForParallel = 0

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
#for parallel
#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior#None
#numericalFluxType = Stokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #None
#numericalFluxType = StokesP_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #None
#auxiliaryVariables=[VelocityAverage()]
auxiliaryVariables=[KTensor()]
#auxiliaryVariables=[BoundaryPressure()]
