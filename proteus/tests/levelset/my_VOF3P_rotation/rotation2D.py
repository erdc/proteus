from proteus import Domain
from proteus.mprans import VOF3P
from proteus import Norms
from proteus import Profiling
import numpy as np

timeIntegration_vof = "SSP33"
#timeIntegration_vof = "FE"

fullNewton=True
# ENTROPY VISCOSITY and ART COMPRESSION PARAMETERS
EDGE_VISCOSITY=1
ENTROPY_VISCOSITY=1
POWER_SMOOTHNESS_INDICATOR=2
LUMPED_MASS_MATRIX=0
FCT=1
cK=0.25
# FOR EDGE BASED ENTROPY VISCOSITY
cE=1.0
cMax=0.1
# FOR SUPG 
shockCapturingFactor_vof=0.2
#Other time parameters
if timeIntegration_vof == "SSP33":
    timeOrder = 3
else:
    timeOrder = 1

runCFL = 0.1#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
lag_shockCapturing_vof=True
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
cDegree_vof=0
pDegree_vof=1
useHex=False
useMetrics=0.0
#
#spatial quadrature orders
rotation_quad_order = 2*pDegree_vof+1
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
#spatial mesh
lRefinement=3
#tag simulation name to level of refinement
nn=nnx=nny=(2**lRefinement)*10+1
nnz=1
he=1.0/(nnx-1.0)
L=[1.0,1.0]

unstructured=False #True for tetgen, false for tet or hex from rectangular grid
box=Domain.RectangularDomain(L=(1.0,1.0),
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
#end time of simulation
T=1.0
#number of output time steps
nDTout = 10
#smoothing factors
#eps
epsFactHeaviside=epsFactDirac=epsFact_vof=1.5 #1.5
#
if useMetrics:
    shockCapturingFactor_vof=0.5
    lag_shockCapturing_vof=True

#use absolute tolerances on al models
atolRedistance = max(1.0e-12,0.1*he)
atolConservation = max(1.0e-12,0.001*he**2)
atolVolumeOfFluid= max(1.0e-12,0.001*he**2)
atolLevelSet     = max(1.0e-12,0.001*he**2)
#controls
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
#redist solver
fmmFlag=0
#
if useHex:
    hex=True
    soname="rotation_c0q"+`pDegree_vof`+"_"+timeIntegration_vof+"_level_"+`lRefinement`
else:
    soname="rotation_c0p"+`pDegree_vof`+"_"+timeIntegration_vof+"_level_"+`lRefinement`

#My Own Coefficients
class MyCoefficients(VOF3P.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[self.modelIndex]
        self.u_dof_old = np.copy(self.model.u[0].dof)
        self.u_dof_old_old = np.copy(self.model.u[0].dof)
        self.velx_tn_dof = np.zeros(self.model.u[0].dof.shape,'d')+1E10
        self.vely_tn_dof = np.zeros(self.model.u[0].dof.shape,'d')+1E10
        self.velz_tn_dof = np.zeros(self.model.u[0].dof.shape,'d')+1E10
        self.flux_plus_dLij_times_soln = np.zeros(self.model.u[0].dof.shape,'d')
        self.q_v = np.zeros((self.model.mesh.nElements_global,self.model.nQuadraturePoints_element,self.model.nSpace_global),'d')+1E10
        self.ebqe_v = np.zeros((self.model.mesh.nExteriorElementBoundaries_global,self.model.nElementBoundaryQuadraturePoints_elementBoundary,self.model.nSpace_global),'d')
        self.model.q[('velocity',0)]=self.q_v
        self.model.ebqe[('velocity',0)]=self.ebqe_v
        self.ebqe_phi = np.zeros(self.model.ebqe[('u',0)].shape,'d') #NOTE: this is not needed (is for LS)
        # Divergence. Assume the velocity is div free
        self.q_div_velocity = np.zeros(self.model.q[('u', 0)].shape,'d')
        self.ebqe_div_velocity = np.zeros(self.model.ebqe[('u', 0)].shape,'d')
        # FOR POROSITY
        self.q_vos = np.zeros(self.model.q[('u', 0)].shape, 'd')
            
    def preStep(self,t,firstStep=False):
        # SAVE OLD SOLUTIONS
        self.u_dof_old_old = np.copy(self.u_dof_old)
        self.u_dof_old = np.copy(self.model.u[0].dof)

        # COMPUTE VELOCITY
        pi = np.pi

        # GET VELOCITY AT DOFs (FOR EDGE BASED METHODS)
        x_dof = self.model.u[0].femSpace.mesh.nodeArray[:,0]
        y_dof = self.model.u[0].femSpace.mesh.nodeArray[:,1]
        self.velx_tn_dof = -2.0*pi*y_dof
        self.vely_tn_dof = 2.0*pi*x_dof
        self.velz_tn_dof = 0.0*x_dof

        # GET VELOCITY AT QUADRATURE POINTS (FOR CELL BASE METHODS)
        x = self.model.q['x'][...,0]
        y = self.model.q['x'][...,1]
        x_boundary = self.model.ebqe['x'][...,0]
        y_boundary = self.model.ebqe['x'][...,1]

        ############
        # ROTATION #
        ############
        self.q_v[...,0]  = -2.0*pi*(y-0.5)
        self.q_v[...,1]  =  2.0*pi*(x-0.5)

        self.ebqe_v[...,0]  = -2.0*pi*(y_boundary-0.5)
        self.ebqe_v[...,1]  =  2.0*pi*(x_boundary-0.5)
        
        #DIVERGENCE OF VEOCITY
        self.q_div_velocity = 0*x
        self.ebqe_div_velocity[...] = 0*x_boundary

        ###################
        # PERIODIC VORTEX #
        ###################
        #T=8
        #yy=y
        #xx=x
        #self.q_v[...,0] = -2*np.sin(pi*yy)*np.cos(pi*yy)*np.sin(pi*xx)**2*np.cos(pi*t/T)
        #self.q_v[...,1] = 2*np.sin(pi*xx)*np.cos(pi*xx)*np.sin(pi*yy)**2*np.cos(pi*t/T)
        
        #self.ebqe_v[...,0]  = -2*np.sin(pi*y_boundary)*np.cos(pi*y_boundary)*np.sin(pi*x_boundary)**2*np.cos(pi*t/T)
        #self.ebqe_v[...,1]  =  2*np.sin(pi*x_boundary)*np.cos(pi*x_boundary)*np.sin(pi*y_boundary)**2*np.cos(pi*t/T)
        
        #DIVERGENCE OF VEOCITY
        #self.q_div_velocity = 0*x
        #self.ebqe_div_velocity[...] = 0*x_boundary
        
        ###############
        # TRANSLATION #
        ###############
        #self.q_v[...,0]  = 0.0
        #self.q_v[...,1]  = -1.0

        # CHECK MASS
        if self.checkMass:
            self.m_pre = Norms.scalarDomainIntegral(self.model.q['dV_last'],
                                                    self.model.q[('m',0)],
                                                    self.model.mesh.nElements_owned)
            Profiling.logEvent("Phase  0 mass before VOF step = %12.5e" % (self.m_pre,),level=2)

        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        if(self.FCT==1):
            self.model.FCTStep()
        if self.checkMass:
            self.m_post = Norms.scalarDomainIntegral(self.model.q['dV'],
                                                     self.model.q[('m',0)],
                                                     self.model.mesh.nElements_owned)
            Profiling.logEvent("Phase  0 mass after VOF step = %12.5e" % (self.m_post,),level=2)
            #Compare mass before and after step
            #np.testing.assert_almost_equal(self.m_pre,self.m_post, err_msg="Mass before and after step are not almost equal", verbose=True)

        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        pass
