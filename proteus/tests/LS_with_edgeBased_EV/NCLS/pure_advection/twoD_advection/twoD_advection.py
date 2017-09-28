from proteus import Domain
from proteus.mprans import NCLS
from proteus import Norms
from proteus import Profiling 
import numpy as np
import math 

SSPOrder = 2
lRefinement=3

problem = 2 #0: dist function, 1: rotation of smoothed disk, 2: rotation of zalesak disk, 3: periodic vortex
epsCoupez=3 # Relative to h. This is used to saturate distance function
LUMPED_MASS_MATRIX=False
STABILIZATION_TYPE=1 #0:SUPG, 1:EV, 2:Based on smoothness indicator 
ENTROPY_TYPE=2 #1: eta=0.5*u^2, 2: logarithmic 
cE=0.1
runCFL = 0.5

if (problem<3):
    T=1.0
    nDTout = 10
else: 
    T=2.0
    nDTout = 20

# SHOCK CAPTURING PARAMETERS
shockCapturingFactor_ncls=0.2
# OTHER TIME PARAMETERS

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
twoD_advection_quad_order = 2*pDegree_ncls+1
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node

#tag simulation name to level of refinement
nn=nnx=(2**lRefinement)*10+1
nny=nnx
nnz=1
he=1.0/(nnx-1.0)

unstructured=False #True for tetgen, false for tet or hex from rectangular grid
box=Domain.RectangularDomain(L=(1.0,1.0),
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
    soname="twoD_advection_level_"+`lRefinement`
else:
    soname="twoD_advection_level_"+`lRefinement`

class MyCoefficients(NCLS.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[0]
        self.q_v = np.zeros(self.model.q[('dH',0,0)].shape,'d')
        self.ebqe_v = np.zeros(self.model.ebqe[('dH',0,0)].shape,'d')
        self.rdModel = self.model

