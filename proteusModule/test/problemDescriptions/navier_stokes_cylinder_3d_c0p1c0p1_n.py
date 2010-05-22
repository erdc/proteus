from pyadh import *
from pyadh.default_n import *
from navier_stokes_cylinder_3d_p import *

# timeIntegration = FLCBDF
# stepController  = FLCBDF_controller
# systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
timeIntegration = BackwardEuler
stepController = Min_dt_controller
systemStepControllerType = SplitOperator.Sequential_MinModelStep
runCFL=0.9
#stepController  = FixedStep
#systemStepControllerType = SplitOperator.Sequential_FixedStep
#timeIntegration = NoIntegration
#stepController  = Newton_controller
rtol_u[1] = 1.0e-2
rtol_u[2] = 1.0e-2
rtol_u[3] = 1.0e-2
atol_u[1] = 1.0e-2
atol_u[2] = 1.0e-2
atol_u[3] = 1.0e-2

tnList = [0.0,0.1*residence_time,T]

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nLevels = 1
triangleOptions="VApq1.25q12ena%21.16e" % (((0.75*inflow_height)**3)/6.0,)#(((0.25*inflow_height)**3)/6.0,)
#triangleOptions="VpAq1.25ena%f" % (((0.05*inflow_height)**3)/6.0,)
#triangleOptions="pAq30Dena%f" % (0.5*(0.25*inflow_height)**2,)
#triangleOptions="pAq30Dena%f" % (0.5*(0.01*inflow_height)**2,)
#triangleOptions="VpAq1.25en"

if useOpt:
    subgridError = NavierStokesASGS_velocity_pressure_opt(coefficients,nd,lag=True,delayLagSteps=1,hFactor=1.0,noPressureStabilization=False)
    shockCapturing = NavierStokes_SC_opt(coefficients,nd,0.25,lag=True)
else:
    subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True)
    shockCapturing = None
#subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd,lag=False)

massLumping = False


#numericalFluxType = None
numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 25

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0#1.0e-6

nl_atol_res = 1.0e-5

matrix = SparseMatrix

#multilevelLinearSolver = LU
multilevelLinearSolver = PETSc

#levelLinearSolver = LU
levelLinearSolver = PETSc

linearSmoother = GaussSeidel

linTolFac = 0.001
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

#conservativeFlux = {0:'pwl',1:'pwl',2:'pwl'}
conservativeFlux = {0:'pwl-bdm'}#,1:'pwl',2:'pwl',3:'pwl'}

#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho),
#                   RecirculationLength(cylinder_center[0]+cylinder_radius)]
#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho),PressureProfile(flag=5,center=cylinder_center,radius=cylinder_radius),RecirculationLength()]
#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho)]
#auxiliaryVariables=[]
