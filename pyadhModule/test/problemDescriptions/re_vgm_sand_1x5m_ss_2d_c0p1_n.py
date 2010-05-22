from pyadh import *
from pyadh.default_n import *
from re_vgm_sand_1x5m_ss_2d_p import *

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
#cek added high order for use with ASGS/SC
elementQuadrature = {'default':SimplexGaussQuadrature(nd,4),
                     ('u',0):SimplexGaussQuadrature(nd,4)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,4),
                             ('u',0):SimplexGaussQuadrature(nd-1,4)}

nn=3
nLevels = 5

#subgridError = None
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=False)

massLumping = False

numericalFluxType = None

shockCapturing = None
#0.1 works, None,0.05,0.5 don't look like they work
shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.1,lag=False) 

multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 1.0e-3

nl_atol_res = 1.0e-8#1.0e-6

maxNonlinearIts = 100#1

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'pwl'}#{0:'sun-gs-rt0'}#{0:'pwc'}#{0:'pwl'}#{0:'point-eval'}#

#if using sun-gs-rt0 need weak dirichlet conditions
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior
