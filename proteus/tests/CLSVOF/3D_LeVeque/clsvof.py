from proteus import Domain
from proteus import Norms
from proteus import Profiling 
from proteus import Context 
from proteus.mprans import CLSVOF
import numpy as np
import math 

ct=Context.Options([
    # General parameters #
    ("T",3.0,"Final time"),
    ("nDTout",30,"Number of time steps to archive"),
    ("refinement",2,"Level of refinement"),
    ("unstructured",False,"Use unstructured mesh. Set to false for periodic BCs"),
    ("useMetrics",0.0,"Compute h(x) via FE interpolation"),
    # parameters for nonlinear problem
    ("timeOrder",2,"Order of time integration"),
    ("lambdaFact",1.0,"parameter for nonlinear regularization. lambda=epsFactDiffusion")
],mutable=True)

cfl=0.33

# number of space dimensions #
nd=3

# General parameters #
parallel = True # if True use PETSc solvers
linearSmoother = None
checkMass = False

# Finite element sapce #
pDegree_clsvof=1
useBernstein=True
useHex=False

# quadrature order #
clsvof_quad_order = 2*pDegree_clsvof+1

# parallel partitioning info #
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node

# create mesh #
nn=nnx=nny=nnz=(2**ct.refinement)*25+1
#nn=nnx=nny=nnz=(2**ct.refinement)*10+1
he=1.0/(nnx-1.0)

unstructured=ct.unstructured #True for tetgen, false for tet or hex from rectangular grid
box=Domain.RectangularDomain(L=(1.0,1.0,1.0),
                             x=(0.0,0.0,0.0),
                             name="box");
box.writePoly("box")
if unstructured:
    domain=Domain.PlanarStraightLineGraphDomain(fileprefix="box")
    domain.boundaryTags = box.boundaryTags
    bt = domain.boundaryTags
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    domain = box

soname="clsvof_level_"+`ct.refinement`

class MyCoefficients(CLSVOF.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[0]
        self.rdModel = self.model
        self.q_v = np.zeros(self.model.q[('grad(u)',0)].shape,'d')
        self.ebqe_v = np.zeros(self.model.ebqe[('grad(u)',0)].shape,'d')
        self.ebqe_phi = np.zeros(self.model.ebqe[('u',0)].shape,'d')#cek hack, we don't need this         
        self.q_v_old = np.zeros(self.model.q[('grad(u)',0)].shape,'d')

