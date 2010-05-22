from pyadh import *
from pyadh.default_n import *
from burgers_riem_2d_p import *

timeOrder =2
nStagesTime = timeOrder

DT=None
runCFL = 0.3 #0.1

limiterType = TimeIntegration.DGlimiterDurlofskyP1Lagrange2d
#{0:TimeIntegration.DGlimiterP1Lagrange2d}
#{0:TimeIntegration.DGlimiterDurlofskyP1Lagrange2d}  #None

timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller
nDTout = 40

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

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
numericalFluxType =  BurgersNumericalFlux#RusanovNumericalFlux_Diagonal#Advection_DiagonalUpwind#BurgersNumericalFlux#

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

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
