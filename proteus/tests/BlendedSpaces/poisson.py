from __future__ import division
from past.utils import old_div
from proteus import Domain
from proteus import Norms
from proteus import Profiling 
from proteus import Context 
from proteus.mprans import BlendedSpaces
import numpy as np
import math 

ct=Context.Options([
    # General parameters #
    ("pDegree",2,"Order of the polynomial approximation"),
    ("refinement",3,"Level of refinement"),
    ("useBernstein",False,"Use Bernstein basis functions"),
    ("useLinearOnQuadraticCube",False,"Use Linear FE space on quadratic mesh")
],mutable=True)

PROBLEM_TYPE=0
AUTOMATED_ALPHA=1
ALPHA_FOR_GALERKIN_SOLUTION=0
# 0: Poisson equation
# 1: Steady transport problem
# 2: Anisotropic diffusion

COMPOSITE_QUAD_RULE=True

# Finite element sapce #
useHex=True

# number of space dimensions #
nd=2

# General parameters #
parallel = False 

# parallel partitioning info #
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node

##########
# DOMAIN #
##########
if PROBLEM_TYPE==0:
    nn=(2**ct.refinement)*10+1
elif PROBLEM_TYPE==1:
    nn=51 #72  ### #111  #105 min that works....  #251 ==> 126
else: 
    nn=21 #19
    nn=21
#
he=1.0/(nn-1.0)

box=Domain.RectangularDomain(L=(1.0,1.0),
                             x=(0.0,0.0),
                             name="box");
box.writePoly("box")
domain = box

soname="poisson"

class MyCoefficients(BlendedSpaces.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[0]
        self.q_v = np.zeros(self.model.q[('grad(u)',0)].shape,'d')
        self.ebqe_v = np.zeros(self.model.ebqe[('grad(u)',0)].shape,'d')
        self.ebqe_phi = np.zeros(self.model.ebqe[('u',0)].shape,'d')#cek hack, we don't need this 
        
