from pyadh import *
from pyadh.default_n import *
from twp_darcy_split_Buckley_Leverett_sat_p import *


limiterType =None
timeOrder = 1
nStagesTime = timeOrder
#type of time integration formula
timeIntegration = BackwardEuler_cfl
#timeIntegration = ForwardEuler_A
#general type of integration (Forward or to SteadyState)
timeIntegrator  = ForwardIntegrator
#stepController = Min_dt_controller
runCFL = 0.6
    
#timeIntegration = LinearSSPRKintegration
#stepController=Min_dt_RKcontroller




DT=None
nDTout = 1 #int(T/DT)
print "nDTout",nDTout

#femSpaces = {0:DG_Constants}
femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=51
nLevels=1
subgridError = None
massLumping=False
shockCapturing = None

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#numericalFluxType = RusanovNumericalFlux_Diagonal_Diffusion_IIPG
numericalFluxType = RusanovNumericalFlux_Diagonal

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 10

fullNewtonFlag = True

tolFac = 0.0#0.001

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None

archiveFlag = ArchiveFlags.EVERY_USER_STEP
