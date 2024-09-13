from proteus import *
from proteus.default_n import *
from md_trans_p import *

timeIntegration = BackwardEuler_cfl
DT = 0.1
runCFL = 0.3
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)


massLumping=False

subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=True)

shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton
fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

matrix = SparseMatrix
maxNonlinearIts = 15


if parallel:                                                                     
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    #have to have a numerical flux in parallel                                   
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior         
    #for true residual test                                                      
    #linearSolverConvergenceTest = 'r-true'                                      
    #to allow multiple models to set different ksp options                       
    linear_solver_options_prefix = 'trans_' 
    linearSmoother = None                                                        
else:                                          
    multilevelLinearSolver = LU                                                  
    levelLinearSolver = LU                                                  
    smoother = Jacobi

linTolFac = 1.0e-10

conservativeFlux = None