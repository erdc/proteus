from pyadh import *
from pyadh.default_n import *
from kEpsilon_step_2d_p import *

#timeIntegration = BackwardEuler_cfl
#stepController  = Min_dt_controller
explicit = False
LDG = True
if explicit:
    runCFL=0.001#0.3
    spaceOrder = 0;   
    timeOrder = min(spaceOrder+1,3)
    nStagesTime=timeOrder
    timeIntegration = SSPRKPIintegration
    stepController=Min_dt_RKcontroller
elif useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    runCFL = 0.1#None
else:
    timeIntegration = FLCBDF
    stepController = FLCBDF_controller_sys
    atol_u[0] = 1.0e-2
    rtol_u[0] = 1.0e-2
    atol_u[1] = 1.0e-2
    rtol_u[1] = 1.0e-2

#runCFL = 0.5
DT = None
nDTout = 1

femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
             1:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,space_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,space_quad_order)

#nn=3
#triangleOptions = "Aq30Dena%f" % (0.15**2 / 6.0)
#triangleOptions = "Aq30Den"
#nLevels = 1


subgridError = None

if LDG:
    numericalFluxType = RusanovLDG
else:
    numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
maxNonlinearIts =10

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
