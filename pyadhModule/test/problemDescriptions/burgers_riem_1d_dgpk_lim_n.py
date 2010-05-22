from pyadh import *
from pyadh.default_n import *
from burgers_riem_1d_p import *

timeOrder =3 # max allowed for Nonlinear SSPRK
nStagesTime = timeOrder

DT=None
runCFL = 0.1
limiterType = TimeIntegration.DGlimiterPkMonomial1d

#mwf BackwardEuler 
timeIntegration =  SSPRKPIintegration
stepController=Min_dt_RKcontroller
nDTout = 10

#recall time integration order must be k+1
#femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
#highest can go with current quadrature
femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}


elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 41
nLevels = 1

subgridError = None

massLumping = False

class BurgersNumericalFlux(ConvexOneSonicPointNumericalFlux):
   def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        ConvexOneSonicPointNumericalFlux.__init__(self,vt,getPointwiseBoundaryConditions,
                                                  getAdvectiveFluxBoundaryConditions,
                                                  getDiffusiveFluxBoundaryConditions,
                                                  sonicPoint=0.0,sonicFlux=0.0)
   #
#

numericalFluxType = RusanovNumericalFlux_Diagonal#BurgersNumericalFlux#RusanovNumericalFlux_Diagonal

shockCapturing = None

multilevelNonlinearSolver  = NLNI

usingSSPRKNewton = True
levelNonlinearSolver = SSPRKNewton#Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
