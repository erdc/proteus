from proteus import *
from proteus.default_n import *
from cylinder_john_p import *

timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep

rtol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
atol_u[1] = 1.0e-2
atol_u[2] = 1.0e-2

#tnList = [0,T]
dt_fixed = 1.0e-2 #2.5e-1
dt_init = 1.0e-2
nDTout = int(T/dt_fixed)
dt_init = min(dt_init,0.5*dt_fixed)
tnList = [0.0,dt_init]+[i*dt_fixed for i in range(1,nDTout+1)]


femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}
# femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              1:C0_AffineQuadraticOnSimplexWithNodalBasis,
#              2:C0_AffineQuadraticOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

# elementQuadrature = SimplexLobattoQuadrature(nd,1)

# elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

nLevels = 1
diameterConstraint = 0.075*inflow_height
areaConstraint = 0.5*diameterConstraint**2
triangleOptions="pAq30Dena%f" % (areaConstraint,)

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True)

numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior

massLumping = False

shockCapturing = None

multilevelNonlinearSolver  = NLNI

levelNonlinearSolver = Newton

maxNonlinearIts = 25

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0#1.0e-6

atol = 1.0e-8

nl_atol_res = 1.0e-8


matrix = SparseMatrix

# multilevelLinearSolver = LU
# levelLinearSolver = LU

multilevelLinearSolver = PETSc
levelLinearSolver = PETSc
parallelPartitioningType = MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 1

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'pwl',1:'pwl',2:'pwl'}

auxiliaryVariables=[BoundaryForce(D=2.0*cylinder_radius,Ubar=Ubar,rho=rho),
                    PressureProfile(flag=5,center=cylinder_center,radius=cylinder_radius),
                    RecirculationLength(cylinder_center[0]+cylinder_radius)]

