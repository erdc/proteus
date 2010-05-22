from pyadh import *
from pyadh.default_n import *
from navier_stokes_duct_3d_p import *

timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
rtol_u[3] = 1.0e-3
atol_u[1] = 1.0e-3
atol_u[2] = 1.0e-3
atol_u[3] = 1.0e-3

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)

#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nnx=81
nny=2
nnz=41
nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True)
#subgridError = StokesASGS_velocity(coefficients,nd)

shockCapturing = None

#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton #no nested iteration for psitc

levelNonlinearSolver = Newton

fullNewtonFlag = True
maxNonlinearIts =20

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc

levelLinearSolver = LU
#levelLinearSolver = PETSc

linTolFac = 0.001

conservativeFlux = None#{0:'pwl'}

