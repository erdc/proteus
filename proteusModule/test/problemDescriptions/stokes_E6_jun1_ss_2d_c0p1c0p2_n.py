from pyadh import *
from pyadh.default_n import *
from stokes_E6_jun1_ss_2d_p import *

timeIntegration = NoIntegration
timeIntegrator  = ForwardIntegrator


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineQuadraticOnSimplexWithNodalBasis,
             2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

#make sure each pore throat has at least one element in interior
triangleOptions="cq32Den"#"cq30Den"#"cq30Dena0.0001"#
nn=3
nLevels = 1

subgridError = None

massLumping = False

shockCapturing = None

numericalFluxType = None

maxNonlinearIts = 2
multilevelNonlinearSolver  = Newton#NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0#1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'pwl',1:'point-eval',2:'point-eval'}#None

#if want weak boundary conditions make sure use StokesP coefficients
numericalFluxType = None#Stokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior
