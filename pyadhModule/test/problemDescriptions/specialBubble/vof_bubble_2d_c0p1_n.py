from pyadh import *
from pyadh.default_n import *
from bubble import *
from vof_bubble_2d_p import *

elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,bubble_quad_order)


timeIntegration = BackwardEuler_cfl
stepController=Min_dt_controller
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
numericalFluxType = None
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)#linear
subgridError = Advection_ASGS(coefficients=coefficients,nd=nd,lag=False)
massLumping = False
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

# timeOrder=1
# nStagesTime = timeOrder
# class SSPRKwrap(LinearSSPRKintegration):
#     """
#     wrap SSPRK so default constructor uses
#     the order I want and runCFL without
#     changing VectorTransport
#     """
#     def __init__(self,vt):
#         LinearSSPRKintegration.__init__(self,vt,timeOrder,runCFL)
#         return
# timeIntegration = SSPRKwrap 
# stepController=Min_dt_RKcontroller
# femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
# subgridError = None
# numericalFluxType = Advection_DiagonalUpwind
# shockCapturing = None
# subgridError = None
# massLumping = False

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

atol = 1.0e-8

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
   
linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
