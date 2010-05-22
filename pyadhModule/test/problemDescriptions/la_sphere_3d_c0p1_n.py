from pyadh import *
from pyadh.default_n import *
from la_sphere_3d_p import *


timeOrder = 2
nStagesTime = timeOrder
runCFL = 0.1
class SSPRKwrap(SSPRKintegration):
    """
    wrap SSPRK so default constructor uses
    the order I want and runCFL without
    changing VectorTransport
    """
    def __init__(self,vt):
        SSPRKintegration.__init__(self,vt,timeOrder,runCFL)
        return
    #
nStagesTime = 1
timeIntegration = BackwardEuler
runCFL=0.5
#timeIntegration = ForwardEuler
#timeIntegration = SSPRKwrap 
#stepController=Min_dt_RKcontroller

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=11
nLevels = 1
DT = None #1.0e-3
nDTout = 10

#readFromTriangle =  False # False True
#triangleFileBase = 'rectmesh.3'
#if readFromTriangle:
#    nLevels = 1


#subgridError = None
subgridError = Advection_ASGS(coefficients,nd,stabFlag='1',lag=True)

massLumping = False

numericalFluxType = None

#shockCapturing = None
shockCapturing = ResGrad_SC(coefficients,nd)

multilevelNonlinearSolver  = NLNI #Newton

levelNonlinearSolver = Newton

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

