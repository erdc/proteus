from pyadh import *
from pyadh.default_n import *
from bubble import *
from ls_bubble_2d_p import *


# class SSPRKwrap(SSPRKPIintegration):
#     """
#     wrap SSPRK so default constructor uses
#     the order I want and runCFL without
#     changing VectorTransport
#     """
#     def __init__(self,vt):
#         SSPRKPIintegration.__init__(self,vt,timeOrder,runCFL)
#         return
#     #
# #end wrapper
# timeIntegrator = ForwardIntegrator
# #timeIntegration = ForwardEuler
# if useBackwardEuler:
#     timeIntegration = BackwardEuler
# #timeIntegration = FLCBDF
# #runCFL = 1.5
# #timeIntegration = FLCBDF
# timeOrder = 1
# nStagesTime = timeOrder
# #timeIntegration = BackwardEuler
# #BackwardEuler,SSPRKwrap
# #timeIntegration = SSPRKwrap 
# #stepController=Min_dt_RKcontroller
# #SSPRK
# #timeOrder = 3
# #nStagesTime = timeOrder
# #runCFL = 0.15
# #DT=1.e-5
# #DT=None
# #nDTout = 100
# timeIntegrator = ForwardIntegrator
# timeIntegration = FLCBDF
# stepController = FLCBDF_controller
if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    rtol_u[0] = 1.0e-2
    atol_u[0] = 1.0e-2
#timeIntegration = FLCBDF
#stepController = FLCBDF_controller_sys
#rtol_u[0] = 1.0e-2
#atol_u[0] = 1.0e-2

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elif spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,bubble_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,bubble_quad_order)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)#H is linear, non need to lag

massLumping = False

numericalFluxType = None
#numericalFluxType = Advection_DiagonalUpwind

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)

multilevelNonlinearSolver  = Newton#NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

atol = 1.0e-8#should be linear with lagging

maxNonlinearIts = 50

matrix = SparseMatrix

if usePETSc:
    numericalFluxType = NF_base # does nothing

    multilevelLinearSolver = PETSc
    
    levelLinearSolver = PETSc
else:
    multilevelLinearSolver = LU
    
    levelLinearSolver = LU
   
linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
