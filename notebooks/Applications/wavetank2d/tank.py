from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent

#  Discretization -- input options

Refinement = 12
genMesh=True
useOldPETSc=False
useSuperlu=False#True
timeDiscretization='be'#'vbdf'#'be','flcbdf'
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 0.0
useOnlyVF = False
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega
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
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
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

# Domain and mesh
L = (0.584,0.350)
#L = (0.584 , 0.584)
obst_portions = (0.024,0.048) #(width,height)
obst_x_start = 0.292 # start x coordinate of the obstacle; caution to be in the domain's range
obst_x_end = obst_x_start + obst_portions[0] # end x coordinate of the obstacle; caution to be in the domain's range
obst = (obst_x_start,obst_portions[1],obst_x_end) #coordinates of the obstacle to be used to define the boundary


he = L[0]/float(4*Refinement-1)
#he*=0.5
#he*=0.5
#he*=0.5
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

structured=False
#structured=True # Trying out a structured mesh

class PointGauges(AV_base):
    def  __init__(self,gaugeLocations={'pressure_1':(0.5,0.5,0.0)}):
        AV_base.__init__(self)
        self.locations=gaugeLocations
        self.flags={}
        self.files={}#will be opened  later
        pointFlag=100
        for name,point in self.locations.iteritems():
            self.flags[name] = pointFlag
            pointFlag += 1
    def attachModel(self,model,ar):
        self.model=model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.p = model.levelModelList[-1].u[0].dof
        self.u = model.levelModelList[-1].u[1].dof
        self.v = model.levelModelList[-1].u[2].dof
        return self
    def attachAuxiliaryVariables(self,avDict):
        return self
    def calculate(self):
        import numpy as  np
        for name,flag  in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name+'.txt','w')
                self.files[name].write('%22.16e %22.16e %22.16e  %22.16e  %22.16e\n' % (self.vertices[vnMask,0],self.vertices[vnMask,1],self.p[vnMask],self.u[vnMask],self.v[vnMask]))

class LineGauges(AV_base):
    def  __init__(self,gaugeEndpoints={'pressure_1':((0.5,0.5,0.0),(0.5,1.8,0.0))},linePoints=10):
        import numpy as  np
        AV_base.__init__(self)
        self.endpoints=gaugeEndpoints
        self.flags={}
        self.linepoints={}
        self.files={}#while open later
        pointFlag=1000
        for name,(pStart,pEnd) in self.endpoints.iteritems():
            self.flags[name] = pointFlag
            p0 = np.array(pStart)
            direction = np.array(pEnd) - p0
            self.linepoints[name]=[]
            for scale in np.linspace(0.0,1.0,linePoints):
                self.linepoints[name].append(p0 + scale*direction)
            pointFlag += 1
    def attachModel(self,model,ar):
        self.model=model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.p = model.levelModelList[-1].u[0].dof
        self.u = model.levelModelList[-1].u[1].dof
        self.v = model.levelModelList[-1].u[2].dof
        return self
    def attachAuxiliaryVariables(self,avDict):
        return self
    def calculate(self):
        import numpy as  np
        for name,flag  in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name+'.txt','w')
                for x,y,p,u,v in zip(self.vertices[vnMask,0],self.vertices[vnMask,1],self.p[vnMask],self.u[vnMask],self.v[vnMask]):
                    self.files[name].write('%22.16e %22.16e %22.16e  %22.16e  %22.16e\n' % (x,y,p,u,v))

class LineGauges_phi(AV_base):
    def  __init__(self,gaugeEndpoints={'pressure_1':((0.5,0.5,0.0),(0.5,1.8,0.0))},linePoints=10):
        import numpy as  np
        AV_base.__init__(self)
        self.endpoints=gaugeEndpoints
        self.flags={}
        self.linepoints={}
        self.files={}#while open later
        pointFlag=1000
        for name,(pStart,pEnd) in self.endpoints.iteritems():
            self.flags[name] = pointFlag
            p0 = np.array(pStart)
            direction = np.array(pEnd) - p0
            self.linepoints[name]=[]
            for scale in np.linspace(0.0,1.0,linePoints):
                self.linepoints[name].append(p0 + scale*direction)
            pointFlag += 1
    def attachModel(self,model,ar):
        self.model=model
        self.vertexFlags = model.levelModelList[-1].mesh.nodeMaterialTypes
        self.vertices = model.levelModelList[-1].mesh.nodeArray
        self.phi = model.levelModelList[-1].u[0].dof
        return self
    def attachAuxiliaryVariables(self,avDict):
        return self
    def calculate(self):
        import numpy as  np
        for name,flag  in self.flags.iteritems():
            vnMask = self.vertexFlags == flag
            if vnMask.any():
                if not self.files.has_key(name):
                    self.files[name] = open(name+'_phi.txt','w')
                for x,y,phi in zip(self.vertices[vnMask,0],self.vertices[vnMask,1],self.phi[vnMask]):
                    self.files[name].write('%22.16e %22.16e %22.16e\n' % (x,y,phi))

pointGauges = PointGauges(gaugeLocations={'pointGauge_pressure':(0.293,0.05,0.0)})
lineGauges  = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825':((0.4,0.0,0.0),(0.4,0.58,0.0))},linePoints=20)
#'lineGauge_x/H=1.653':((0.99,0.0,0.0),(0.99,1.8,0.0))
lineGauges_phi  = LineGauges_phi(lineGauges.endpoints,linePoints=20)


if useHex:
    nnx=4*Refinement+1
    nny=2*Refinement+1
    hex=True
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=4*Refinement
        nny=2*Refinement
    else:
        vertices=[[0.0,0.0],#0
                  [obst[0],0.0], #1
                  [obst[0],obst[1]], #2
                  [obst[2],obst[1]], #3
                  [obst[2],0.0],#4
                  [L[0],0.0],#5
                  [L[0],L[1]],#6
                  [0.0,L[1]]] #7


        vertexFlags=[boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['bottom'],
                     boundaryTags['top'],
                     boundaryTags['top']]



        segments=[[0,1],
                  [1,2],
                  [2,3],
                  [3,4],
                  [4,5],
                  [5,6],
                  [6,7],
                  [7,0]]


        segmentFlags=[boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['bottom'],
                      boundaryTags['right'],
                      boundaryTags['top'],
                      boundaryTags['left']]



        regions=[[0.1*L[0] ,0.5*L[1]]]
        regionFlags=[1]
        #for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
        #    vertices.append(gaugeCoordinates)
        #    vertexFlags.append(pointGauges.flags[gaugeName])
        #for gaugeName,gaugeLines in lineGauges.linepoints.iteritems():
        #    for gaugeCoordinates in gaugeLines:
        #        vertices.append(gaugeCoordinates)
        #        vertexFlags.append(lineGauges.flags[gaugeName])
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions="VApq30Dena%8.8f" % ((he**2)/2.0,)

logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
# Time stepping
T=0.3
dt_fixed = 0.02
dt_init = min(0.1*dt_fixed,0.001)
runCFL=0.9
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = False#True
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
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.1
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.1
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
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
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-8,0.01*he**2)
vof_nl_atol_res = max(1.0e-8,0.01*he**2)
ls_nl_atol_res = max(1.0e-8,0.01*he**2)
rd_nl_atol_res = max(1.0e-8,0.01*he)
mcorr_nl_atol_res = max(1.0e-8,0.01*he**2)
kappa_nl_atol_res = max(1.0e-8,0.01*he**2)
dissipation_nl_atol_res = max(1.0e-8,0.01*he**2)

#turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,-9.8]

# Initial condition
waterLine_x =0.146
waterLine_z = 0.292

def signedDistance(x):
    phi_x = x[0]-waterLine_x
    phi_z = x[1]-waterLine_z
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)
