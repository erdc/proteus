from pyadh import *
from pyadh.default_n import *
from impes_sat_p import *

#type of time integration formula
#timeIntegration = BackwardEuler
timeIntegration = FLCBDF
stepController = FLCBDF_controller_sys
rtol_u[0] = 1.0e-4
atol_u[0] = 1.0e-4
#timeOrder=1
#nTimeStages=1
#timeIntegration = LinearSSPRKIntegration
#stepController = Min_dt_RKcontroller
#runCFL = 0.33
#limiterType =TimeIntegration.DGlimiterP1Lagrange1d
    

runCFL=0.1

#nDTout=100

DT=1.0e1 
nDTout = int(T/DT)
print "nDTout",nDTout

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

#femSpaces = {0:DG_Constants}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=101
nLevels=1
subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
massLumping=False
shockCapturing = None

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
numericalFluxType = RusanovNumericalFlux_Diagonal_Diffusion_IIPG
#numericalFluxType = Advection_DiagonalUpwind

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 100

fullNewtonFlag = True

tolFac = 0.0#0.001

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None

archiveFlag = ArchiveFlags.EVERY_USER_STEP
