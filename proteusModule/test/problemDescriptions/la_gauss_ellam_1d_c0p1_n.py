from pyadh import *
from pyadh.default_n import *
from la_gauss_ellam_1d_p import *

#BackwardEuler
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
DT = 1.0e-3
runCFL = 4.5#0.2
nDTout = 5
dtout = 0.5*T/nDTout
tmpt=numpy.arange(nDTout+1)*dtout
tnList=list(tmpt)
tmpt += 0.5*T
tnList.extend(list(tmpt[1:]))


femSpaces = dict((i,C0_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

quad_order = 5#5
elementQuadrature = SimplexGaussQuadrature(nd,quad_order)#CompositeTrapezoidalEdge(quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order) #1 point regardless
    
nn = 41#101
nLevels = 1
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
shockCapturing = None

analyticalTracking = 2
if velocityIsTransient:
    if analyticalTracking and velocitySpaceFlag == 'rt0':
        particleTracking = Tracking.LinearAdvection_RT0Velocity_PT123A
    else:
        particleTracking = Tracking.LinearAdvection_C0P1Velocity_PT123
else:
    if analyticalTracking and velocitySpaceFlag == 'rt0':
        particleTracking = Tracking.SteadyState_LinearAdvection_RT0Velocity_PT123A
    else:
        particleTracking = Tracking.SteadyState_LinearAdvection_C0P1Velocity_PT123
    #particleTracking = Tracking.SteadyState_LinearAdvection_C0P1Velocity_PT123
    #mwf try RT0 
    #particleTracking = Tracking.SteadyState_LinearAdvection_RT0Velocity_PT123A
    
particleTracking_params = {}
for ci in range(nc):
    particleTracking_params[('atol_tracking',ci)]=1.0e-7
    particleTracking_params[('rtol_tracking',ci)]=0.0
    particleTracking_params[('sf_tracking',ci)]=0.9
    particleTracking_params[('dn_safe_tracking',ci)]=1.0e-7
#how to evaluate old mass integral
useBackwardTrackingForOldMass = False#True

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

