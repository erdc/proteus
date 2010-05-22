from pyadh import *
from pyadh.default_n import *
#make sure to change name of imported p file when you copy this file
from navier_stokes_poiseuille_ss_2d_p import *

timeIntegration = NoIntegration

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,5)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

nn=51
#number of refinement  levels, might want to decrease for 3D

nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts =100

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU
levelLinearSolver = LU
#multilevelLinearSolver = PETSc
#levelLinearSolver = PETSc

linTolFac = 0.001

conservativeFlux = None
#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior

#conservativeFlux = {0:'pwl',1:'pwl',2:'pwl'}
auxiliaryVariables=[VelocityAverage()]
