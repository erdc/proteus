from pyadh import *
from pyadh.default_n import *
from navier_stokes_cavity_2d_p import *

#timeIntegration = FLCBDF
#stepController  = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
#timeIntegration = BackwardEuler_cfl
#stepController  =  Min_dt_controller
#systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
#restrictFineSolutionToAllMeshes = True
rtol_u[1] = 1.0e-2
rtol_u[2] = 1.0e-2
atol_u[1] = 1.0e-2
atol_u[2] = 1.0e-2
#tnList = [0.0,0.1*residence_time,T]
tnList = [0.0,T]

timeIntegration = NoIntegration
stepController  = Newton_controller

timeIntegration = PsiTCtte_new
#tnList = [0.0,0.1*L[0]/speed]
stepController  = PsiTCtte_controller
rtol_res[0] = 0.0
rtol_res[1] = 0.0
rtol_res[2] = 0.0
atol_res[0] = 1.0e-8
atol_res[1] = 1.0e-8
atol_res[2] = 1.0e-8
#atol_u[0] = 1.0e-3
#rtol_u[0] = 1.0e-3
#timeIntegration = BackwardEuler
#stepController = FixedStep
#runCFL = 0.9

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

# femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

# elementQuadrature = SimplexGaussQuadrature(nd,5)

# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

nn=321
nn=11
nLevels = 1
#DT=1.0e2
#nDTout = int(T/DT)
#DT = None
nDTout = 1
#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False,delayLagSteps=2)
subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False,delayLagSteps=0)
#subgridError = StokesASGS_velocity(coefficients,nd)

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior #need weak for parallel and global conservation

shockCapturing = None

#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton #no nested iteration for psitc

levelNonlinearSolver = Newton

fullNewtonFlag = True

maxNonlinearIts =100
maxLineSearches=0
tolFac = 0.0

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = PETSc#LU

levelLinearSolver = PETSc# LU

linTolFac = 0.001

conservativeFlux = None#{0:'pwl'}

