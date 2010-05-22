from pyadh import *
from pyadh.default_n import *
from sloshbox3d import *
from ls_consrv_sloshbox_3d_p import *


timeIntegrator = ForwardIntegrator
timeIntegration = BackwardEuler#NoIntegration
timeIntegration = NoIntegration

if spaceOrder == 1:
    femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elif spaceOrder == 2:
    femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
#femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}
#femSpaces = {0:DG_Constants}
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,sloshbox_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,sloshbox_quad_order)

# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

subgridError = None
#subgridError = HamiltonJacobi_ASGS(coefficients,nd,lag=False)

massLumping = False

numericalFluxType = DoNothing#Diffusion_IIPG_exterior#None
#numericalFluxType = Advection_Diagonal_average
#numericalFluxType = Advection_DiagonalUpwind
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
if usePETSc:
    numericalFluxType = DoNothing#Diffusion_IIPG_exterior
shockCapturing = None

#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-4#1.0e-8#*(he**3/6.0)#1.0e-4

maxNonlinearIts = 50

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

if usePETSc:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    
linearSmoother = GaussSeidel

linTolFac = 1.0e-6

conservativeFlux = None
auxiliaryVariables = [AuxiliaryVariables.ConservationHistoryMC("sloshbox3d"+`lRefinement`+"p"+`spaceOrder`)]
