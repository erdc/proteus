from pyadh import *
from pyadh.default_n import *
from re_vgm_sand_10x10m_ss_2d_p import *

timeIntegration = NoIntegration
class SteadyStateIntergratorWrap(SteadyStateIntegrator):
    def __init__(self,mlvtran,mlnl,dtMeth,nOptions,ssatol=1.0e-3,ssrtol=1.0e-3,
                 maxsteps=500,stepExact=True,ignoreFailure=False):
        SteadyStateIntegrator.__init__(self,mlvtran,mlnl,dtMeth,nOptions,ssatol,ssrtol,
                                       maxsteps,stepExact,ignoreFailure)
    #def
#class
#timeIntegration = PsiTCtte
#timeIntegrator = SteadyStateIntergratorWrap

DT = 1.0

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#appears to work with NLNI or Newton but not subgridErro
#elementQuadrature = {'default':SimplexLobattoQuadrature(nd,3),
#                     ('u',0):SimplexGaussQuadrature(nd,3)}
#elementBoundaryQuadrature ={'default':SimplexGaussQuadrature(nd-1,2),
#                            ('u',0):SimplexGaussQuadrature(nd-1,2)}
#doesn't work for testStuff velPP v2
#elementBoundaryQuadrature = {'default':SimplexLobattoQuadrature(nd-1,3),
#                             ('u',0):SimplexGaussQuadrature(nd-1,3)}

#cek added high order for use with ASGS/SC
elementQuadrature = {'default':SimplexGaussQuadrature(nd,3),
                     ('u',0):SimplexGaussQuadrature(nd,3)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3),
                             ('u',0):SimplexGaussQuadrature(nd-1,3)}

nn=3
nLevels = 4

subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)
#to try with PTC
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = None
#for satform
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)
#for headform
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)
#try with PTC
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)
#shockCapturing = ResGradQuad_SC(coefficients,nd,shockCapturingFactor=2.0,lag=False)
#shockCapturing = ResGradJuanes_SC(coefficients,nd,shockCapturingFactor=4.0,uSC=4.5,lag=False)

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

maxNonlinearIts = 100#1

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

#conservativeFlux = {0:'point-eval'}
conservativeFlux = {0:'pwl'}
#conservativeFlux = None

