"""
Bar floating in half-filled tank
"""
from math import *
from proteus import *
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.MeshAdaptPUMI import MeshAdaptPUMI

Profiling.verbose=True
from proteus import Context
opts=Context.Options([
    ("bar_dim", (0.33,0.33,0.2), "Dimensions of the bar"),
    ("tank_dim", (1.0,1.0,1.0), "Dimensions of the tank"),
    ("water_surface_height",0.5,"Height of free surface above bottom"),
    ("speed",0.0,"Speed of current if non-zero"),
    ("bar_height",0.85,"Initial height of bar center above bottom"),
    ("bar_rotation",(0,0,0),"Initial rotation about x,y,z axes"),
    ("refinement_level",0,"Set maximum element diameter to he/2**refinement_level"),
    ("gen_mesh",False,"Generate new mesh"),
    ("T",0.1,"Simulation time"),
    ("dt_init",0.001,"Initial time step"),
    ("cfl",0.33,"Target cfl"),
    ("nsave",5,"Number of time steps to  save"),
    ("parallel",True,"Run in parallel"),
    ("free_x",(0.0,0.0,1.0),"Free translations"),
    ("free_r",(1.0,1.0,0.0),"Free rotations")])

#----------------------------------------------------
# Physical properties
#----------------------------------------------------
rho_0=998.2
nu_0 =1.004e-6

rho_1=1.205
nu_1 =1.500e-5

sigma_01=0.0

g=[0.0,0.0,-9.81]

#----------------------------------------------------
# Domain - mesh - quadrature
#----------------------------------------------------
nd = 3

(bar_length,bar_width,bar_height)  = opts.bar_dim

L=opts.tank_dim

x_ll = (0.0,0.0,0.0)

waterLevel   =  opts.water_surface_height

bar_center = (0.5*L[0],0.5*L[1],opts.bar_height)
speed = opts.speed

#set up barycenters for force calculation
barycenters = numpy.zeros((8,3),'d')
barycenters[7,:] = bar_center

bar_mass    = bar_length*bar_width*bar_height*0.5*(rho_0+rho_1)

bar_cg      = [0.0,0.0,0.0]

bar_inertia = [[(L[1]**2+L[2]**2)/12.0, 0.0                    , 0.0                   ],
               [0.0                   , (L[0]**2+L[2]**2)/12.0 , 0.0                   ],
               [0.0                   , 0.0                    , (L[0]**2+L[1]**2)/12.0]]

RBR_linCons  = [1,1,0]
RBR_angCons  = [1,0,1]


nLevels = 1

he = 0.1 #coarse grid
he *=(0.5)**opts.refinement_level

boundaries = [ 'bottom', 'front', 'right', 'back', 'left', 'top', 'obstacle' ]
boundaryTags = dict([(key,i+1) for (i,key) in enumerate(boundaries)])

#Floating
faceMap = { 'front'     : [82],
            'right'     : [24],
            'back'      : [78],
            'left'      : [42],
            'bottom'    : [80],
            'top'       : [76],
            'obstacle'  : [163,178,158,173,153,168] }

#Falling
faceMap = { 'front'     : [82],
            'right'     : [24],
            'back'      : [78],
            'left'      : [42],
            'bottom'    : [80],
            'top'       : [76],
            'obstacle'  : [174,164,179,154,159,169] }

faceList = []
for boundary in boundaries:
  if boundary in faceMap:
    faceList.append(faceMap[boundary])
  else:
    faceList.append([])

domain = Domain.PUMIDomain(dim=3)
domain.faceList = faceList

nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element

adaptMesh = True
adaptMesh_nSteps = 1
adaptMesh_numIter = 3

domain.PUMIMesh = MeshAdaptPUMI.MeshAdaptPUMI(hmax=0.1,hmin=0.05,
                                              numIter=adaptMesh_numIter,sfConfig="ERM",maType="isotropic")
comm = Comm.init()

#Remember to change bar height: Floating = 0.25, Falling=0.85
#Remember to change Face Maps @_@
case_dir = "simModel/falling_bar"
model_dir = "%s-Procs" % comm.size()
case_mesh = "Floating_Bar_Falling_coarse.smb"
case_model= "Floating_Bar_Falling_coarse.smd"
input_model= "%s/%s" % (case_dir,case_model)
input_mesh = "%s/%s/%s" % (case_dir,model_dir,case_mesh)
domain.PUMIMesh.loadModelAndMesh(input_model, input_mesh)

restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 0

quad_order = 3

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
openTop = False
openSides = False
openEnd = True
smoothBottom = False
smoothObstacle = False
movingDomain=True
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init=opts.dt_init
T = opts.T
nDTout=opts.nsave
dt_out =  (T-dt_init)/nDTout
runCFL = opts.cfl

#----------------------------------------------------
water_depth  = waterLevel-x_ll[2]

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not opts.parallel
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

#  Discretization
nd = 3
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,3)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
	basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
	basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


# Numerical parameters
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.5
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.5
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-12,0.001*he**2)
vof_nl_atol_res = max(1.0e-12,0.001*he**2)
ls_nl_atol_res = max(1.0e-12,0.001*he**2)
mcorr_nl_atol_res = max(1.0e-12,0.0001*he**2)
rd_nl_atol_res = max(1.0e-12,0.01*he)
kappa_nl_atol_res = max(1.0e-12,0.001*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.001*he**2)
mesh_nl_atol_res = max(1.0e-12,0.001*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x,t):
    p_L = 0.0
    phi_L = L[2] - waterLevel
    phi = x[2] - waterLevel
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

import ode

def near_callback(args, geom1, geom2):
    """Callback function for the collide() method.

    This function checks if the given geoms do collide and
    creates contact joints if they do.
    """

    # Check if the objects do collide
    contacts = ode.collide(geom1, geom2)

    # Create contact joints
    world,contactgroup = args
    for c in contacts:
        c.setBounce(0.2)
        c.setMu(5000)
        j = ode.ContactJoint(world, contactgroup, c)
        j.attach(geom1.getBody(), geom2.getBody())

class RigidBar(AuxiliaryVariables.AV_base):
    def __init__(self,density=1.0,bar_center=(0.0,0.0,0.0),bar_dim=(1.0,1.0,1.0),barycenters=None,he=1.0,cfl_target=0.9,dt_init=0.001):
        self.dt_init = dt_init
        self.he=he
        self.cfl_target=cfl_target
        self.world = ode.World()
        #self.world.setERP(0.8)
        #self.world.setCFM(1E-5)
        self.world.setGravity(g)
        self.g = np.array([g[0],g[1],g[2]])
        self.space = ode.Space()
        eps_x = L[0]- 0.75*L[0]
        eps_y = L[1]- 0.75*L[1]
        #tank geometry
        #self.tankWalls = [ode.GeomPlane(self.space, (1,0,0) ,x_ll[0]+eps_x),
        #                  ode.GeomPlane(self.space, (-1,0,0),-(x_ll[0]+L[0]-eps_x)),
        #ode.GeomPlane(self.space, (0,1,0) ,x_ll[1]+eps_y),
        #                  ode.GeomPlane(self.space, (0,-1,0) ,-(x_ll[1]+L[1]-eps_y))]
        #mass/intertial tensor of rigid bar
        #eps_x = L[0]- 0.45*L[0]
        #eps_y = L[1]- 0.45*L[1]
        #self.tankWalls = [ode.GeomPlane(space, (1,0,0) ,x_ll[0]+eps_x),
        #                  ode.GeomPlane(space, (-1,0,0),-(x_ll[0]+L[0]-eps_x)),
        #                  ode.GeomPlane(space, (0,1,0) ,x_ll[1]+eps_y),
        #                  ode.GeomPlane(space, (0,-1,0) ,-(x_ll[1]+L[1]-eps_y))]
        #self.tank = ode.GeomBox(self.space,(0.45,0.45,0.3))
        #self.tank.setPosition((0.5,0.5,0.6))
        #self.contactgroup = ode.JointGroup()
        self.M = ode.Mass()
        self.totalMass = density*bar_dim[0]*bar_dim[1]*bar_dim[2]
        self.M.setBox(density,bar_dim[0],bar_dim[1],bar_dim[2])
        #bar body
        self.body = ode.Body(self.world)
        self.body.setMass(self.M)
        self.body.setFiniteRotationMode(1)
        #bar geometry
        self.bar = ode.GeomBox(self.space,bar_dim)
        self.bar.setBody(self.body)
        self.bar.setPosition(bar_center)
        self.boxsize = (bar_dim[0],bar_dim[1],bar_dim[2])
        #contact joints
        self.contactgroup = ode.JointGroup()
        self.last_position=bar_center
        self.position=bar_center
        self.last_velocity=(0.0,0.0,0.0)
        self.velocity=(0.0,0.0,0.0)
        self.h=(0.0,0.0,0.0)
        self.rotation = np.eye(3)
        self.last_rotation = np.eye(3)
        self.last_rotation_inv = np.eye(3)
        self.barycenters=barycenters
        self.init=True
        self.bar_dim = bar_dim
        self.last_F = np.zeros(3,'d')
        self.last_M = np.zeros(3,'d')

    def attachModel(self,model,ar):
        self.model=model
        self.ar=ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        assert(flagMin >= 0)
        assert(flagMax <= 7)
        self.nForces=flagMax+1
        assert(self.nForces <= 8)
        return self
    def get_u(self):
        return self.last_velocity[0]
    def get_v(self):
        return self.last_velocity[1]
    def get_w(self):
        return self.last_velocity[2]
    def calculate_init(self):
        self.last_F = None
        self.calculate()
    def calculate(self):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = self.dt_init
        F = self.model.levelModelList[-1].coefficients.netForces_p[7,:] + self.model.levelModelList[-1].coefficients.netForces_v[7,:];
        M = self.model.levelModelList[-1].coefficients.netMoments[7,:]
        logEvent("x Force " +`self.model.stepController.t_model_last`+" "+`F[0]`)
        logEvent("y Force " +`self.model.stepController.t_model_last`+" "+`F[1]`)
        logEvent("z Force " +`self.model.stepController.t_model_last`+" "+`F[2]`)
        logEvent("x Moment " +`self.model.stepController.t_model_last`+" "+`M[0]`)
        logEvent("y Moment " +`self.model.stepController.t_model_last`+" "+`M[1]`)
        logEvent("z Moment " +`self.model.stepController.t_model_last`+" "+`M[2]`)
        logEvent("dt " +`dt`)
        scriptMotion=False
        linearOnly=False
        if self.last_F == None:
            self.last_F = F.copy()
        if scriptMotion:
            velocity = np.array((0.0,0.3/1.0,0.0))
            logEvent("script pos="+`(np.array(self.position)+velocity*dt).tolist()`)
            self.body.setPosition((np.array(self.position)+velocity*dt).tolist())
            self.body.setLinearVel(velocity)
        else:
            if linearOnly:
                Fstar = 0.5*(F+self.last_F) + np.array(self.world.getGravity())
                velocity_last = np.array(self.velocity)
                velocity = velocity_last + Fstar*(dt/self.totalMass)
                velocity[0] = 0.0
                vmax = self.he*self.cfl_target/dt
                vnorm = np.linalg.norm(velocity,ord=2)
                if vnorm > vmax:
                    velocity *= vmax/vnorm
                    logEvent("Warning: limiting rigid body velocity from "+`vnorm`+" to "+`vmax`)
                position_last = np.array(self.position)
                position = position_last + 0.5*(velocity_last + velocity)*dt
                self.body.setPosition(position.tolist())
                self.body.setLinearVel(velocity.tolist())
                msg = """
Fstar         = {0}
F             = {1}
F_last        = {2}
dt            = {3:f}
velocity      = {4}
velocity_last = {5}
position      = {6}
position_last = {7}""".format(Fstar,F,self.last_F,dt,velocity,velocity_last,position,position_last)
                logEvent(msg)
            else:
                nSteps=10
                # vnorm = np.linalg.norm((F+self.g)*dt)
                # vmax = self.he*self.cfl_target/dt
                # if vnorm > vmax:
                #     F *= vmax/vnorm
                #     logEvent("Warning: limiting rigid body velocity from "+`vnorm`+" to "+`vmax`)

                Fstar=F#0.5*(F+self.last_F)
                Mstar=M#0.5*(M+self.last_M)
                for i in range(nSteps):
                    self.body.setForce((Fstar[0]*opts.free_x[0],
                                        Fstar[1]*opts.free_x[1],
                                        Fstar[2]*opts.free_x[2]))
                    self.body.setTorque((Mstar[0]*opts.free_r[0],
                                         Mstar[1]*opts.free_r[1],
                                         Mstar[2]*opts.free_r[2]))
                    #self.space.collide((self.world,self.contactgroup), near_callback)
                    self.world.step(dt/float(nSteps))
        #self.contactgroup.empty()
        self.last_F[:] = F
        self.last_M[:] = M
        x,y,z = self.body.getPosition()
        u,v,w = self.body.getLinearVel()
        self.barycenters[7,0]=x
        self.barycenters[7,1]=y
        self.barycenters[7,2]=z
        self.last_velocity=copy.deepcopy(self.velocity)
        self.last_position=copy.deepcopy(self.position)
        self.last_rotation=self.rotation.copy()
        self.last_rotation_inv = inv(self.last_rotation)
        self.position=(x,y,z)
        self.velocity=(u,v,w)
        self.rotation=np.array(self.body.getRotation()).reshape(3,3)
        self.h = (self.position[0]-self.last_position[0],
                  self.position[1]-self.last_position[1],
                  self.position[2]-self.last_position[2])
        logEvent("%1.2fsec: pos=(%21.16e, %21.16e, %21.16e) vel=(%21.16e, %21.16e, %21.16e) h=(%21.16e, %21.16e, %21.16e)" % (self.model.stepController.t_model_last,
                                                                                    self.position[0],
                                                                                    self.position[1],
                                                                                    self.position[2],
                                                                                    self.velocity[0],
                                                                                    self.velocity[1],
                                                                                                            self.velocity[2],
                                                                                                            self.h[0],
                                                                                                            self.h[1],
                                                                                                            self.h[2]))

bar = RigidBar(density=0.5*(rho_0+rho_1),bar_center=bar_center,bar_dim=opts.bar_dim,barycenters=barycenters,he=he,cfl_target=0.9*opts.cfl,dt_init=opts.dt_init)
