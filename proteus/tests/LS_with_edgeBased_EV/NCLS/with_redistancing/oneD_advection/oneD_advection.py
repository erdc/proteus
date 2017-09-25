from proteus import Domain
from proteus.mprans import NCLS
from proteus import Norms
from proteus import Profiling 
import numpy as np
import math 

orderSSP = 2
lRefinement=3
T=1.0
nDTout = 10
runCFL = 0.5

problem=0 #0: distance function, 1: saturated distance function
pure_redistancing=True
COUPEZ = True
DO_SMOOTHING = False
DO_REDISTANCING = True
STABILIZATION_TYPE=1 #Just for transport equation. 1: EV, 2: smoothness based

# COUPEZ and REDISTANCING
redist_tolerance=0.1
epsCoupez=3 #For saturation
epsFactRedistance=0.33 #For signed dist function
lambda_coupez = 0.1 #Strength of redistancing and coupez force

if problem==0: #dist function
    SATURATED_LEVEL_SET = False
    ENTROPY_TYPE=1 # eta(u)=0.5*u**2
    cE=1.0
else: #saturated distance function
    SATURATED_LEVEL_SET = True
    ENTROPY_TYPE=2 
    cE=0.1

# SHOCK CAPTURING PARAMETERS
shockCapturingFactor_ncls=0.2

lag_shockCapturing_ncls=True
#if True uses PETSc solvers
parallel = False
linearSmoother = None
#compute mass balance statistics or not
checkMass=False
#number of space dimensions
nd=2
#time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0:1.0e-4}
atol_u = {0:1.0e-4}
rtol_res = {0:1.0e-4}
atol_res = {0:1.0e-4}
#
#spatial approximation orders
cDegree_ncls=0 
pDegree_ncls=1
useBernstein=False
useHex=False
useMetrics=0.0
#
#spatial quadrature orders
oneD_advection_quad_order = 2*pDegree_ncls+1
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node

#tag simulation name to level of refinement
nn=nnx=(2**lRefinement)*10+1
nny=(nnx-1)/10+1
nnz=1
he=1.0/(nnx-1.0)

unstructured=False #True for tetgen, false for tet or hex from rectangular grid
box=Domain.RectangularDomain(L=(1.0,0.1),
                             x=(0.0,0.0),
                             name="box");
box.writePoly("box")
if unstructured:
    domain=Domain.PlanarStraightLineGraphDomain(fileprefix="box")
    domain.boundaryTags = box.boundaryTags
    bt = domain.boundaryTags
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    domain = box

#eps
epsFactHeaviside=epsFactDirac=1.5
if useMetrics:
    shockCapturingFactor=0.5
    lag_shockCapturing=True

#use absolute tolerances on al models
atolLevelSet     = max(1.0e-12,0.001*he**2)
#controls
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
#redist solver
fmmFlag=0
#
if useHex:
    hex=True

soname="oneD_advection_level_"+`lRefinement`

class MyCoefficients(NCLS.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[0]
        self.q_v = np.zeros(self.model.q[('dH',0,0)].shape,'d')
        self.ebqe_v = np.zeros(self.model.ebqe[('dH',0,0)].shape,'d')
        self.rdModel = self.model
        
