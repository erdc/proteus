from proteus import Domain
from proteus.mprans import NCLS
from proteus import Norms
from proteus import Profiling 
import numpy as np
import math 

#timeIntegration_ncls = "SSP33" 
timeIntegration_ncls = "FE" 
lRefinement=3
T=0.5
nDTout = 10

pure_redistancing=False
redist_tolerance=0.1
epsCoupez=2
epsFactRedistance=0.33
lambda_coupez = 1.0
# ENTROPY VISCOSITY 
EDGE_VISCOSITY=1
ENTROPY_VISCOSITY=1
LUMPED_MASS_MATRIX=0
# SHOCK CAPTURING PARAMETERS
shockCapturingFactor_ncls=0.2
# OTHER TIME PARAMETERS
if timeIntegration_ncls == "SSP33":
    timeOrder = 3
else:
    timeOrder = 1

runCFL = 0.25
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
#soname="oneD_advectioncgp2_bdf2_mc"+`lRefinement`
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
    from tank2dDomain import *
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
    soname="oneD_advection_c0q"+`pDegree_ncls`+"_"+timeIntegration_ncls+"_"+`timeOrder`+"_level_"+`lRefinement`
else:
    soname="oneD_advection_c0p"+`pDegree_ncls`+"_"+timeIntegration_ncls+"_"+`timeOrder`+"_level_"+`lRefinement`

class MyCoefficients(NCLS.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[0]
	self.u_old_dof = np.copy(self.model.u[0].dof)
	self.u_old_dof_old = np.copy(self.model.u[0].dof)
        self.q_v = np.zeros(self.model.q[('dH',0,0)].shape,'d')+1E10
        self.ebqe_v = np.zeros(self.model.ebqe[('dH',0,0)].shape,'d')
        self.model.q[('velocity',0)]=self.q_v
        self.model.ebqe[('velocity',0)]=self.ebqe_v
        if self.RD_modelIndex != None:
            #print self.RD_modelIndex,len(modelList)
            self.rdModel = modelList[self.RD_modelIndex]
        else:
            self.rdModel = self.model

    def preStep(self,t,firstStep=False):
        # SAVE OLD SOLUTIONS
        self.u_old_dof_old = np.copy(self.u_old_dof)
        self.u_old_dof = np.copy(self.model.u[0].dof)

        self.q_v[...,0] = 1
        self.q_v[...,1] = 0 
        self.ebqe_v[...,0] = 1
        self.ebqe_v[...,1] = 0

        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        pass
