from pyadh import *
from pyadh.default_n import *
from impes_sat_p import *

#general type of integration (Forward or to SteadyState)
timeIntegrator = ForwardIntegrator
timeIntegration = FLCBDF
stepController = FLCBDF_controller_sys
rtol_u[1] = 1.0e-4
atol_u[1] = 1.0e-4
runCFL=None#0.001
# timeOrder =1
# nStagesTime = timeOrder
# class SSPRKwrap(SSPRKintegration):
#     """
#     wrap SSPRK so default constructor uses
#     the order I want and runCFL without
#     changing VectorTransport
#     """
#     def __init__(self,vt):
#         SSPRKintegration.__init__(self,vt,timeOrder,runCFL)
#         return
#     #
# #end wrapper
    
# #mwf BackwardEuler 
# timeIntegration = SSPRKwrap 
# stepController=Min_dt_RKcontroller

#nDTout=100

DT=1.0e1 
nDTout = int(T/DT)
print "nDTout",nDTout
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

#elementQuadrature = SimplexLobattoQuadrature(nd,1)
#elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)

#nn=3
#nLevels = 1
nn=101
nLevels=1
subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)


massLumping=False
#massLumping=True

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.9,lag=True)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.75,lag=False)
#shockCapturing = ScalarAdvection_SC(coefficients,nd,shockCapturingFactor=0.25,lag=False)

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = Newton
maxNonlinearIts = 25

fullNewtonFlag = True

tolFac = 0.0#0.001

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None
