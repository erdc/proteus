from pyadh import *
from pyadh.default_n import *
from navier_stokes_cylinder_2d_p import *
timeIntegrator = ForwardIntegrator
useFLCBDF = False#True#False
useGustafsson = True#False#True#

nDTout = 100#1001#2000#int(T/DT)
if useFLCBDF:
    timeIntegration = FLCBDF
    stepController  = FLCBDF_controller
    systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
    rtol_u[1] = 1.0e-2
    rtol_u[2] = 1.0e-2
    atol_u[1] = 1.0e-2
    atol_u[2] = 1.0e-2
elif useGustafsson:
    #timeIntegration = BackwardEuler
    timeIntegration = VBDF
    timeOrder = 2
    stepController = Min_dt_cfl_controller#GustafssonFullNewton_dt_controller#FixedStep
    #nDTout = 1#int(T/DT)#int(T/DT) #100#int(T/DT)
    #for controlling time stepping
    nonlinearConvergenceRate_ref = 0.2#0.3#0.2
    useInitialGuessPredictor= True
    stepExact = True
    systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
    runCFL = 0.1
    rtol_u[1] = 1.0e-2#1.0e-2
    rtol_u[2] = 1.0e-2#1.0e-2
    atol_u[1] = 1.0e-2#1.0e-2
    atol_u[2] = 1.0e-2#1.0e-2
else:
    timeIntegration = BackwardEuler_cfl
    stepController  = Min_dt_controller
    runCFL = 0.1#None
    #timeIntegration = BackwardEuler
    
    #DT=1.0e-1
    #DT=0.0000001
    #stepController = HeuristicNL_dt_controller#FixedStep
    #nDTout = 1#int(T/DT)#int(T/DT) #100#int(T/DT)
    #for controlling time stepping
    #nonlinearIterationsFloor = 4
    #nonlinearIterationsCeil  = 8
    #dtNLgrowFactor = 2
    #dtNLreduceFactor = 0.5
    #dtNLfailureReduceFactor = 0.5
    #useInitialGuessPredictor= True
    #stepExact = True
    if T/float(nDTout) > 1.0e-2:
        tnList = [0.0,1.0e-2]
    else:
        tnList = [0.0]
    tnList.extend([n*float(T/nDTout) for n in range(1,nDTout+1)])

#tnList = [0.0,0.1*residence_time,T]

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}
elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
if useOpt:
    subgridError = NavierStokesASGS_velocity_pressure_opt(coefficients,nd,lag=True,delayLagSteps=1,hFactor=1.0,noPressureStabilization=False)
    shockCapturing = NavierStokes_SC_opt(coefficients,nd,0.25,lag=True)
else:
    subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=0,hFactor=1.0)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)


nn=3
nLevels = 1
triangleOptions="pAq30Dena%f" % (0.5*(0.1*inflow_height)**2,)
#triangleOptions="pAq30Dena%f" % (0.5*(0.05*inflow_height)**2,)
#triangleOptions="pAq30Dena%f" % (0.5*(0.05*inflow_height)**2,)
#triangleOptions="pAq30Dena%f" % (0.5*(0.0125*inflow_height)**2,)
#triangleOptions="pAq30Dena%f" % (0.5*(0.0025*inflow_height)**2,)
#triangleoptions="Aq30Den"

#triangleOptions="pAq30Dena%f" % (0.5*(0.025*inflow_height)**2,)
#triangleOptions="Aq30Den"

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True,noPressureStabilization=False)
#subgridError = NavierStokesTransientSubScalesASGS_velocity_pressure(coefficients,nd,lag=True,delayLagSteps=0,noPressureStabilization=True,
#                                                                    trackSubScales=True,
#                                                                    limit_tau_t=True,tau_t_limit_max=0.9,tau_t_limit_min=0.1)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)
#numericalFluxType = None
numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation
#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior #need weak for parallel and global conservation
#numericalFluxType = Stokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

massLumping = False

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 5#25

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0#1.0e-6

nl_atol_res = 1.0e-6

matrix = SparseMatrix
parallel = True
if parallel:
    #for petsc do things lie
    #"-ksp_type bcgsl -pc_type asm -pc_asm_type basic -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
    #-pc_type lu -pc_factor_mat_solver_package
    #-ksp_type bcgsl -pc_type asm -pc_asm_type basic -sub_ksp_type preonly -sub_pc_type lu -sub_pc_factor_mat_solver_package superlu_dist -ksp_rtol 0.0 -ksp_atol 1.0e-10
    #can also set -pc_asm_overlap 2 with default asm type (restrict)
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    #pick number of layers to use in overlap 
    nLayersOfOverlapForParallel = 1
    #type of partition
    parallelPartitioningType = MeshParallelPartitioningTypes.node
    #parallelPartitioningType = MeshParallelPartitioningTypes.element
    #numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior
    numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_SIPG_exterior
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

#conservativeFlux = {0:'pwl',1:'pwl',2:'pwl'}
#conservativeFlux = {0:'point-eval',1:'point-eval',2:'point-eval'}
#conservativeFlux = {0:'pwl',1:'pwl',2:'pwl'}

#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho),
#                    PressureProfile(flag=5,center=cylinder_center,radius=cylinder_radius),
#                    RecirculationLength(cylinder_center[0]+cylinder_radius)]
#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho),PressureProfile(flag=5,center=cylinder_center,radius=cylinder_radius)]
#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho)]
#auxiliaryVariables=[rc]
#auxiliaryVariables=[AverageVelocity()]
#auxiliaryVariables=[BoundaryPressure()]
#parallelPartitioningType = MeshParallelPartitioningTypes.node
