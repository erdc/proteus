from pyadh import *
from pyadh.default_n import *
from la_gauss_ellam_3d_p import *

#BackwardEuler
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
runCFL = 8.5#8.5#0.2
nDTout = 5
dtout = 0.5*T/nDTout
tmpt=numpy.arange(nDTout+1)*dtout
tnList=list(tmpt)
tmpt += 0.5*T
tnList.extend(list(tmpt[1:]))

femSpaces = dict((i,C0_AffineLinearOnSimplexWithNodalBasis) for i in range(nc))

quad_order = 2
#elementQuadrature =CompositeTrapezoidalTriangle(quad_order)
elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order) #1 point regardless

if modelColumn:
    nnx = 2; nny = 2; nnz = 21
else:
    nn = 21#11
nLevels = 1
numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
shockCapturing = None

#tracking options
zeroSolutionTol_track = {0:1.0e-7}
particleTracking_params = {}
if velocitySpaceFlag == 'rt0':
    analyticalTracking = 0
    if analyticalTracking == 1:
        particleTracking = Tracking.SteadyState_LinearAdvection_RT0Velocity_AnalyticalTracking_2d
        for ci in range(nc):
            particleTracking_params[('zeroTol',ci)]=1.0e-6
            particleTracking_params[('localVelocityRepresentationFlag',ci)]=2
    elif analyticalTracking == 2:
        if velocityIsTransient:
            particleTracking = Tracking.LinearAdvection_RT0Velocity_PT123A
        else:
            particleTracking = Tracking.SteadyState_LinearAdvection_RT0Velocity_PT123A
        for ci in range(nc):
            particleTracking_params[('atol_tracking',ci)]=1.0e-7
            particleTracking_params[('rtol_tracking',ci)]=0.0
            particleTracking_params[('sf_tracking',ci)]=0.9
            particleTracking_params[('dn_safe_tracking',ci)]=1.0e-7
            particleTracking_params[('localVelocityRepresentationFlag',ci)]=2
    else:
        if  velocityIsTransient:
            particleTracking = Tracking.LinearAdvection_RT0Velocity_PT123
        else:
            particleTracking = Tracking.SteadyState_LinearAdvection_RT0Velocity_PT123
        for ci in range(nc):
            particleTracking_params[('atol_tracking',ci)]=1.0e-7
            particleTracking_params[('rtol_tracking',ci)]=0.0
            particleTracking_params[('sf_tracking',ci)]=0.9
            particleTracking_params[('dn_safe_tracking',ci)]=1.0e-7
            particleTracking_params[('localVelocityRepresentationFlag',ci)]=2
elif velocitySpaceFlag == 'bdm1':
    if velocityIsTransient:
        particleTracking = Tracking.LinearAdvection_BDM1Velocity_PT123
    else:
        particleTracking = Tracking.SteadyState_LinearAdvection_BDM1Velocity_PT123
    for ci in range(nc):
        particleTracking_params[('atol_tracking',ci)]=1.0e-7
        particleTracking_params[('rtol_tracking',ci)]=0.0
        particleTracking_params[('sf_tracking',ci)]=0.9
        particleTracking_params[('dn_safe_tracking',ci)]=1.0e-7
else:
    if velocityIsTransient:
        particleTracking = Tracking.LinearAdvection_C0P1Velocity_PT123
    else:
        particleTracking = Tracking.SteadyState_LinearAdvection_C0P1Velocity_PT123
    for ci in range(nc):
        particleTracking_params[('atol_tracking',ci)]=1.0e-7
        particleTracking_params[('rtol_tracking',ci)]=0.0
        particleTracking_params[('sf_tracking',ci)]=0.9
        particleTracking_params[('dn_safe_tracking',ci)]=1.0e-7
#how to evaluate old mass integral
useBackwardTrackingForOldMass = False
massLumping = False#True

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton
#maxNonlinearIts = 1

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

