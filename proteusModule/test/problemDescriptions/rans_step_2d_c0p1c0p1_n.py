from pyadh import *
from pyadh.default_n import *
from rans_step_2d_p import *

if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    if timeOrder == 2:
        timeIntegration = VBDF
        #stepController = Min_dt_cfl_controller
        stepController = GustafssonFullNewton_dt_controller
        rtol_u[1] = 1.0e-2
        rtol_u[2] = 1.0e-2
        atol_u[1] = 1.0e-2#1.0e-3
        atol_u[2] = 1.0e-2#1.0e-3
        
#    timeIntegration = FLCBDF
#    stepController = FLCBDF_controller_sys
#    rtol_u[1] = 1.0e-2
#    rtol_u[2] = 1.0e-2
#    atol_u[1] = 1.0e-2#1.0e-3
#    atol_u[2] = 1.0e-2#1.0e-3
else:
    timeIntegration = FLCBDF
    stepController  = FLCBDF_controller_sys
    rtol_u[1] = 1.0e-4
    rtol_u[2] = 1.0e-4
    atol_u[1] = 1.0e-4
    atol_u[2] = 1.0e-4


runCFL = 0.1#None
#DT=0.0000001
nDTout = 10#int(T/DT)
DT = T/nDTout #None#T/100.0
tnList = [0,0.01]
tnList.extend([i*DT for i  in range(1,nDTout+1)])

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,space_quad_order)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,space_quad_order)

#nn=3
#0.15**2 / downstream_length
#triangleOptions = "Aq30Dena%f" % (0.15**2 / 6.0)#(0.15**2 / 6.0)
#triangleOptions = "Aq30Den"
#nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=lag_ns_subgridError)
#NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False) #steady-state Re > 100
#NavierStokesASGS_velocity_pressure(coefficients,nd,lag=True) #transient
#subgridError = StokesASGS_velocity_pressure(coefficients,nd)

massLumping = False

shockCapturing = NavierStokes_SC(coefficients,nd,ns_shockCapturingFactor,lag=True)

#numericalFluxType = NavierStokes_Advection_DiagonalUpwind_Diffusion_IIPG_exterior #need weak for parallel and global conservation
numericalFluxType = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

maxNonlinearIts = 20#100#1000

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = atol_ns#1.0e-6

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = {0:'pwl'}
