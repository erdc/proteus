"""
Basic setup for VOF3P tests following Alistair's template
"""
#selling out
from proteus import iproteus
from proteus.iproteus import *
from proteus import default_p as p
from proteus import default_n as n
from proteus import default_s,default_so
from proteus.mprans import VOF3P
import numpy as np
import proteus as pr
reload(p)
reload(n)

import os,sys
sys.path.insert(0,os.pardir)

import rotation2D as rot2D

"""
p file
"""
p.name=rot2D.soname+"_vof"
p.nd = 2
p.L = [1.0,1.0]
p.domain = Domain.RectangularDomain(L=p.L,
                                    x=(0.0,0.0),
                                    name="box");

p.T = rot2D.T
p.LevelModelType = VOF3P.LevelModel

p.coefficients = rot2D.MyCoefficients(epsFact=rot2D.epsFactHeaviside,checkMass=rot2D.checkMass,useMetrics=rot2D.useMetrics,ME_model=0,
                                      EDGE_VISCOSITY=rot2D.EDGE_VISCOSITY, 
                                      ENTROPY_VISCOSITY=rot2D.ENTROPY_VISCOSITY,
                                      cK=rot2D.cK,cE=rot2D.cE,cMax=rot2D.cMax)
                                        
from proteus.ctransportCoefficients import smoothedHeaviside
class init_cond:
    def __init__(self,center=[0.5,0.75,0.5],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1];
        dBubble = self.radius - np.sqrt(dx**2 + dy**2)
        return smoothedHeaviside(rot2D.epsFactHeaviside*he,dBubble)
    def uOfXT(self,X,t):
        return self.uOfX(X)

p.analyticalSolutions = None

def getDBC(x,flag):
    #return lambda x,t: 1.0
    pass

p.dirichletConditions = {0:getDBC}
p.initialConditions  = {0:init_cond(center=[0.5,0.75],radius=0.15)}
p.fluxBoundaryConditions = {0:'outFlow'}

#cek made no flux since v.n = 0 for this v
def getAFBC(x,flag):
   return lambda x,t: 0.0

p.advectiveFluxBoundaryConditions =  {0:getAFBC}
p.diffusiveFluxBoundaryConditions = {0:{}}


"""
n file
"""
n.multilevelNonlinearSolver  = NonlinearSolvers.NLNI
n.levelNonlinearSolver = NonlinearSolvers.Newton
n.fullNewtonFlag = rot2D.fullNewton
n.updateJacobian = False


n.nDTout = rot2D.nDTout
n.DT = p.T/float(n.nDTout)
n.timeIntegration = VOF3P.RKEV
n.stepController = StepControl.Min_dt_RKcontroller
if rot2D.timeIntegration_vof == "SSP33":
    n.timeOrder = 3
    n.nStagesTime = 3
else:
    n.timeOrder = 1
    n.nStagesTime = 1
 
n.runCFL = rot2D.runCFL
n.rtol_u = rot2D.rtol_u
n.atol_u = rot2D.atol_u
n.rtol_res = rot2D.rtol_res
n.atol_res = rot2D.atol_res

#spatial mesh
lRefinement=3
#tag simulation name to level of refinement
n.nn=n.nnx=n.nny=(2**lRefinement)*10+1
n.nnz=1
he=1.0/(n.nnx-1.0)

if rot2D.pDegree_vof==2:
    n.femSpaces = {0:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis}
else:
    n.femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis}

n.subgridError = None #Advection_ASGS(coefficients,nd,lag=False)
n.shockCapturing = VOF3P.ShockCapturing(p.coefficients,p.nd,shockCapturingFactor=rot2D.shockCapturingFactor_vof,
                                      lag=rot2D.lag_shockCapturing_vof)

if rot2D.parallel or p.LevelModelType == VOF3P.LevelModel:
    n.numericalFluxType = NumericalFlux.Advection_DiagonalUpwind_IIPG_exterior

n.elementQuadrature = Quadrature.SimplexGaussQuadrature(p.nd,rot2D.rotation_quad_order)
n.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(p.nd-1,rot2D.rotation_quad_order)

#parallel partitioning info
n.partitioningType = MeshTools.MeshParallelPartitioningTypes.node

n.tolFac = 0.0
n.linTolFac = n.tolFac
n.linearSolverConvergenceTest = 'r-true'

n.nl_atol_res = 100*rot2D.atolVolumeOfFluid
n.l_atol_res = rot2D.atolVolumeOfFluid

n.maxNonlinearIts = 20
n.maxLineSearches = 0

n.matrix = LinearAlgebraTools.SparseMatrix

if rot2D.parallel:
    n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py#PETSc
    n.levelLinearSolver = LinearSolvers.KSP_petsc4py#PETSc
    n.linear_solver_options_prefix = 'vof_'
    n.linearSolverConvergenceTest = 'r-true'
else:
    n.multilevelLinearSolver = LinearSolvers.LU
    n.levelLinearSolver = LinearSolvers.LU

n.conservativeFlux = {}

n.auxiliaryVariables = [AuxiliaryVariables.MassOverRegion()]

so = default_so
so.name = rot2D.soname
so.sList=[default_s]
so.systemStepExact = True

so.needEBQ_GLOBAL  = False
so.needEBQ = False
so.systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
so.archiveFlag = Archiver.ArchiveFlags.EVERY_USER_STEP
so.DT = n.DT
so.tnList = [i*n.DT for i  in range(n.nDTout+1)]
so.useOneArchive = True

"""

"""
opts = None
#ns = NumericalSolution.NS_base(so,[p],[n],so.sList,iproteus.opts)
