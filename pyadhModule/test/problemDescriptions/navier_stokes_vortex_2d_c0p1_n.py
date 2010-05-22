from pyadh import *
from pyadh.default_n import *
from navier_stokes_vortex_2d_p import *

useFLCBDF = True
if useFLCBDF:
    timeIntegration = FLCBDF
    stepController  = FLCBDF_controller
    systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
    restrictFineSolutionToAllMeshes = False
    #rtol_u[1] = 1.0e-1
    #rtol_u[2] = 1.0e-1
    #atol_u[1] = 1.0e-1
    #atol_u[2] = 1.0e-1
    #nl_atol_res = 1.0e-6
    #rtol_u[1] = 1.0e-2
    #rtol_u[2] = 1.0e-2
    #atol_u[1] = 1.0e-2
    #atol_u[2] = 1.0e-2
    #nl_atol_res = 1.0e-4
    #rtol_u[1] = 1.0e-3
    #rtol_u[2] = 1.0e-3
    #atol_u[1] = 1.0e-3
    #atol_u[2] = 1.0e-3
    #nl_atol_res = 1.0e-5
    # rtol_u[1] = 1.0e-4
    # rtol_u[2] = 1.0e-4
    # atol_u[1] = 1.0e-4
    # atol_u[2] = 1.0e-4
    # nl_atol_res = 1.0e-6
    rtol_u[1] = 1.0e-5
    rtol_u[2] = 1.0e-5
    atol_u[1] = 1.0e-5
    atol_u[2] = 1.0e-5
    nl_atol_res = 1.0e-7
    # rtol_u[1] = 1.0e-6
    # rtol_u[2] = 1.0e-6
    # atol_u[1] = 1.0e-6
    # atol_u[2] = 1.0e-6
    # nl_atol_res = 1.0e-8
    tnList = [0.0,T]
else:
    timeIntegration = BackwardEuler
    stepController  = FixedStep
    nDTout = 10
    tnList = [i*T/float(nDTout) for i in range(nDTout+1)]

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

# femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

# elementQuadrature = SimplexGaussQuadrature(nd,5)

# elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

nn=41
nLevels = 1

#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=5)
#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=0)
subgridError = NavierStokesTransientSubScalesASGS_velocity_pressure(coefficients,nd,lag=False,delayLagSteps=0,noPressureStabilization=False,
                                                                    trackSubScales=True,
                                                                    limit_tau_t=False,tau_t_limit_max=0.9,tau_t_limit_min=0.1)

#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

#shockCapturing = ResGradQuadDelayLag_SC(coefficients,nd,shockCapturingFactor=0.05,lag=True,nStepsToDelay=5)

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

maxNonlinearIts =100
maxLineSearches=10
tolFac = 0.0

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = None#{0:'pwl'}

