from pyadh import *
from pyadh.default_n import *
from navier_stokes_cylinder_ss_2d_p import *

#timeIntegration = FLCBDF
#stepController  = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
# timeIntegration = BackwardEuler
# stepController  = FixedStep
#systemStepControllerType = SplitOperator.Sequential_FixedStep
timeIntegration = NoIntegration
stepController  = Newton_controller
rtol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
atol_u[1] = 1.0e-3
atol_u[2] = 1.0e-3



#runCFL = 0.1#None
#DT=1.0e-1
#DT=0.0000001
#nDTout = 100#2000#int(T/DT)

#DT = T/10.0
tnList = [0.0,T]
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False,hFactor=1.0)

nLevels = 1
triangleOptions="pAq30Dena%f" % (0.5*(0.075*inflow_height)**2,)

numericalFluxType = None
numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation

massLumping = False

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 50

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0#1.0e-6

nl_atol_res = 1.0e-6

linTolFac  = 0.0
l_atol_res = 1.0e-8

matrix = SparseMatrix

#multilevelLinearSolver = LU
#levelLinearSolver = LU
#multilevelLinearSolver = PETSc
#levelLinearSolver = PETSc
#"-ksp_type cg -pc_type asm -pc_asm_type basic -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor_draw" or
#-pc_type lu -pc_factor_mat_solver_package
multilevelLinearSolver = KSP_petsc4py
levelLinearSolver = KSP_petsc4py

linearSmoother = GaussSeidel

linTolFac = 0.001

#conservativeFlux = {0:'pwl',1:'pwl',2:'pwl'}
#conservativeFlux = {0:'point-eval',1:'point-eval',2:'point-eval'}

#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho),
#                    PressureProfile(flag=5,center=cylinder_center,radius=cylinder_radius),
#                    RecirculationLength(cylinder_center[0]+cylinder_radius)]
#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho),PressureProfile(flag=5,center=cylinder_center,radius=cylinder_radius),RecirculationLength()]
#auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho)]
#auxiliaryVariables=[]
