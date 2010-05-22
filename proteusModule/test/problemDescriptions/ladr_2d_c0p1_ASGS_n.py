from pyadh import *
from pyadh.default_n import *
from ladr_2d_p import *

useBackwardEuler = True
if useBackwardEuler:
    timeIntegration = BackwardEuler_cfl
    stepController = Min_dt_controller
    DT = None
    runCFL = 0.1
else:
    #timeIntegration = FLCBDF
    #stepController = FLCBDF_controller
    #systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
    timeIntegration = VBDF
    stepController = GustafssonFullNewton_dt_controller
    rtol_u[0] = 1.0e-2
    atol_u[0] = 1.0e-2
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='1',lag=True)
#need shock capturing for asgs looks like
#shockCapturing=ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.35,lag=True)
subgridError = AdvectionDiffusionReactionTransientSubscales_ASGS(coefficients,nd,stabFlag='1',lag=True,trackSubScales=True,useHarariDirectly=False)
#shockCapturing=ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=True)

nnx=11; nny=11
tnList=[float(i)/10.0 for i in range(11)]
matrix = SparseMatrix
multilevelLinearSolver = LU
