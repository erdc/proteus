"""
A helper module for doing air/water flow around the Wigley hull form comparing to Sakamoto et al. 2007
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


hull_length = 1.0
hull_beam   = hull_length/10.0
hull_draft  = hull_length/16.0

L=(4.0*hull_length,
   2.05*hull_length,
   0.75*hull_length)
x_ll = (-1.5*hull_length,
         -L[1]/2.0,
         0.0)

waterLevel   = 0.5*hull_length

hull_center = (0.0,
               0.0,
               0.5*hull_length)


#debug
# L=(1.5*hull_length,
#  3.0*hull_beam,
#  3.0*hull_draft)

# x_ll = (-0.75*hull_length,
#         -L[1]/2.0,
#         0.0)

# L=(3.0*hull_length,
#    3.05*hull_length,
#    3.0*hull_draft)

# x_ll = (-0.75*hull_length,
#          -L[1]/2.0,
#          0.0)

# waterLevel   = 1.5*hull_draft

# hull_center = (0.0,
#               0.0,
#               waterLevel)

#set up barycenters for force calculation
barycenters = numpy.zeros((8,3),'d')
barycenters[7,:] = hull_center
#cek todo these are the same as the 5414, something's not right
hull_mass    = 532.277
hull_cg      = [2.7618104935392300,  0.0 ,0.27953462008339180  ]
hull_inertia = [[28.2823,  0.0,       20.86855 ],
                [0.0    ,  1126.799,    0.0    ],
		[20.86855,  0.0,      1136.371 ]]

RBR_linCons  = [1,1,0]
RBR_angCons  = [1,0,1]


nLevels = 1

he = hull_draft/0.75 #32
#he *= 0.5 #
#he *= 0.5 #
min_he_poly = he
#he *= 0.5 # tetgen fails
#he *=0.5 #4 way on diamond, 8 way on garnet 256-1024 mpi tasks
#he *=0.5 #2048 - mesh3206851
#he *=0.5 #2048 - mesh3206851
#vessel = 'wigley-gmsh'
#genMesh=False
vessel = 'wigley'
genMesh=True#False
#vessel = 'cube'
#vessel = None
#genMesh=True
boundaryTags = { 'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5, 'top':6, 'obstacle':7}
if vessel is 'wigley-gmsh':
    domain = Domain.MeshTetgenDomain(fileprefix="mesh")
    domain.boundaryTags = boundaryTags
else:
    vertices=[[x_ll[0],x_ll[1],x_ll[2]],#0
              [x_ll[0]+L[0],x_ll[1],x_ll[2]],#1
              [x_ll[0]+L[0],x_ll[1]+L[1],x_ll[2]],#2
              [x_ll[0],x_ll[1]+L[1],x_ll[2]],#3
              [x_ll[0],x_ll[1],x_ll[2]+L[2]],#4
              [x_ll[0]+L[0],x_ll[1],x_ll[2]+L[2]],#5
              [x_ll[0]+L[0],x_ll[1]+L[1],x_ll[2]+L[2]],#6
              [x_ll[0],x_ll[1]+L[1],x_ll[2]+L[2]]]#7
    vertexFlags=[boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left'],
                 boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left']]
    facets=[[[0,1,2,3]],
            [[0,1,5,4]],
            [[1,2,6,5]],
            [[2,3,7,6]],
            [[3,0,4,7]],
            [[4,5,6,7]]]
    facetFlags=[boundaryTags['bottom'],
                boundaryTags['front'],
                boundaryTags['right'],
                boundaryTags['back'],
                boundaryTags['left'],
                boundaryTags['top']]
    regions=[[x_ll[0]+0.5*L[0],x_ll[1]+0.5*L[1],x_ll[2]+0.5*L[2]]]
    regionFlags=[1.0]
    holes=[]
    if vessel is 'wigley':
        from math import log
        he_hull = max(min_he_poly, he)
        #print he,he_hull
        #he_hull = he
        n_points_length = int(ceil(hull_length/he_hull))+1
        n_points_draft  = 2*int(ceil(hull_draft/he_hull))+1
        #print "points",n_points_length,n_points_draft
        dx = hull_length/float(n_points_length-1)
        dz = 2.0*hull_draft/float(n_points_draft-1)
        #print "he",he,dx,dz
        #grid on right half of hull
        for i in range(n_points_length):
            for j in range(n_points_draft):
                x = i*dx - 0.5*hull_length
                z = j*dz - hull_draft
                zStar = min(0.0,z)
                y = 0.5*hull_beam*(1.0 - 4.0*(x/hull_length)**2) * (1.0 - (zStar/hull_draft)**2)
                vertices.append([x+hull_center[0],
                                 y+hull_center[1],
                                 z+hull_center[2]])
                vertexFlags.append(boundaryTags['obstacle'])
        def vN_right(i,j):
            return 8 + i*n_points_draft+j
        for i in range(n_points_length-1):
            for j in range(n_points_draft-1):
                if i < n_points_length/2:
                    facets.append([[vN_right(i,j),vN_right(i+1,j+1),vN_right(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_right(i,j),vN_right(i,j+1),vN_right(i+1,j+1)]])
                    facetFlags.append(boundaryTags['obstacle'])
                else:
                    facets.append([[vN_right(i,j),vN_right(i,j+1),vN_right(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_right(i,j+1),vN_right(i+1,j+1),vN_right(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
        #grid on left half of hull
        for i in range(1,n_points_length-1):
            for j in range(1,n_points_draft):
                x = i*dx - 0.5*hull_length
                z = j*dz - hull_draft
                zStar = min(0.0,z)
                y = 0.5*hull_beam*(1.0 - 4.0*(x/hull_length)**2) * (1.0 - (zStar/hull_draft)**2)
                vertices.append([x+hull_center[0],
                                 hull_center[1] - y,
                                 z+hull_center[2]])
                vertexFlags.append(boundaryTags['obstacle'])
        def vN_left(i,j):
            if i== 0 or j==0:
                return vN_right(i,j)
            if i == (n_points_length-1):# or j==(n_points_draft-1):
                return vN_right(i,j)
            else:
                return 8 + n_points_length*n_points_draft+(i-1)*(n_points_draft-1)+j-1
        for i in range(n_points_length-1):
            for j in range(n_points_draft-1):
                if i < n_points_length/2:
                    facets.append([[vN_left(i,j),vN_left(i+1,j+1),vN_left(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_left(i,j),vN_left(i,j+1),vN_left(i+1,j+1)]])
                    facetFlags.append(boundaryTags['obstacle'])
                else:
                    facets.append([[vN_left(i,j),vN_left(i,j+1),vN_left(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
                    facets.append([[vN_left(i,j+1),vN_left(i+1,j+1),vN_left(i+1,j)]])
                    facetFlags.append(boundaryTags['obstacle'])
        topFacet=[]
        for i in range(n_points_length):
            topFacet.append(vN_right(i,n_points_draft-1))
        for i in range(n_points_length-2,0,-1):
            topFacet.append(vN_left(i,n_points_draft-1))
        facets.append([topFacet])
        facetFlags.append(boundaryTags['obstacle'])
        #for v in vertices: print v
        #for f in facets: print f
        holes.append(hull_center)
    if vessel is 'cube':
        nStart = len(vertices)
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] - 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] - 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] + 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        vertices.append([hull_center[0] + 0.5*hull_length,
                         hull_center[1] - 0.5*hull_beam,
                         hull_center[2] + 0.5*hull_draft])
        vertexFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart,nStart+1,nStart+2,nStart+3]])#1
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart,nStart+1,nStart+5,nStart+4]])#2
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+1,nStart+2,nStart+6,nStart+5]])#3
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+2,nStart+3,nStart+7,nStart+6]])#4
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+3,nStart,nStart+4,nStart+7]])#5
        facetFlags.append(boundaryTags['obstacle'])
        facets.append([[nStart+4,nStart+5,nStart+6,nStart+7]])#6
        facetFlags.append(boundaryTags['obstacle'])
        holes.append(hull_center)
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 vertexFlags=vertexFlags,
                                                 facets=facets,
                                                 facetFlags=facetFlags,
                                                 regions=regions,
                                                 regionFlags=regionFlags,
                                                 holes=holes)
    #go ahead and add a boundary tags member
    domain.boundaryTags = boundaryTags
    if vessel:
        domain.writePoly("mesh_"+vessel)
    else:
        domain.writePoly("meshNoVessel")
    triangleOptions="VApq1.35q12feena%21.16e" % ((he**3)/6.0,)
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0

quad_order = 3

#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
openTop = True
openSides = False
openEnd = True
smoothBottom = False
smoothObstacle = False
movingDomain=False#True
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
Fr = 0.25
#Fr = 0.51
#Fr = 0.0
Um = Fr*sqrt(fabs(g[2])*hull_length)
Re = hull_length*Um/nu_0
weak_bc_penalty_constant = 100.0#Re

if Um > 0.0:
    residence_time = hull_length/Um
else:
    residence_time = 1.0
dt_init=0.001
T = 10*residence_time
nDTout=100
#debut
#T = 0.01#
#nDTout=3
dt_out =  (T-dt_init)/nDTout
runCFL = 0.9

#RANS bc info
kInflow = 0.003*Um*Um
#----------------------------------------------------
# Airy wave functions
#----------------------------------------------------
wave_length = 1.5     * hull_length
wave_height = 0.002   * wave_length
wave_angle  = 0.0     * pi/180.0

#----------------------------------------------------
water_depth  = waterLevel-x_ll[2]
wave_length  = 2.0*pi/wave_length
wave_periode = sqrt(-g[2]*wave_length*tanh(wave_length/water_depth))
wave_vel_amp = wave_periode*(wave_height/sinh(wave_length*water_depth))

xy   = lambda x:   cos(wave_angle)*x[0] + sin(wave_angle)*x[1]
kxwt = lambda x,t: wave_length*(xy(x) - Um*t) - wave_periode*t
kzh  = lambda x:   wave_length*min(x[2]-x_ll[2],water_depth)

#================================================
#  Boundary conditon  lambdas
#================================================
u_wave   = lambda x,t: wave_vel_amp * cosh(kzh(x)) * cos(kxwt(x,t)) * cos(wave_angle)  + Um
v_wave   = lambda x,t: wave_vel_amp * cosh(kzh(x)) * cos(kxwt(x,t)) * sin(wave_angle)
w_wave   = lambda x,t: wave_vel_amp * sinh(kzh(x)) * sin(kxwt(x,t))
noslip   = lambda x,t: 0.0
ls_wave  = lambda x,t: -(wave_height * cos(kxwt(x,t)) + waterLevel - x[2])
vof_wave = lambda x,t: 1.0 if ls_wave(x,t) > 0.0 else 0.0

#================================================
# Print run data
#================================================
logEvent("""
Reynolds number    = %16.21e
Froude number      = %16.21e
Hull Speed[M/S]    = %16.21e
Hull flow time[S]  = %16.21e

Wave length[M]     = %16.21e
Wave height[M]     = %16.21e
Wave angle[Rad]    = %16.21e
Wave periode[Hz]   = %16.21e
Wave velocity[M/S] = %16.21e
T                  = %16.21e
nDTout             = %i
""" % (Re,
       Fr,
       Um,
       residence_time,
       wave_length,
       wave_height,
       wave_angle,
       wave_periode,
       wave_vel_amp,
       T,
       nDTout))

#  Discretization -- input options
useOldPETSc=False
useSuperlu = False # set to False if running in parallel with petsc.options
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 0.0
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
    redist_Newton = True
    kappa_shockCapturingFactor = 0.5
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.5
    dissipation_lag_shockCapturing = True
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
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-12,0.01*he**2)
vof_nl_atol_res = max(1.0e-12,0.01*he**2)
ls_nl_atol_res = max(1.0e-12,0.01*he**2)
mcorr_nl_atol_res = max(1.0e-12,0.01*he**2)
rd_nl_atol_res = max(1.0e-12,0.01*he)
kappa_nl_atol_res = max(1.0e-12,0.01*he**2)
dissipation_nl_atol_res = max(1.0e-12,0.01*he**2)
mesh_nl_atol_res = max(1.0e-12,0.01*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

#wave/current properties
windspeed_u = Um
windspeed_v = 0.0
windspeed_w = 0.0

outflowHeight = waterLevel
outflowVelocity = (Um,0.0,0.0)

inflowHeightMean = waterLevel
inflowVelocityMean = (Um,0.0,0.0)

def waveHeight(x,t):
    return waterLevel

def waveVelocity_u(x,t):
    return Um

def waveVelocity_v(x,t):
    return 0.0

def waveVelocity_w(x,t):
    return 0.0

def wavePhi(x,t):
    return x[2] - waveHeight(x,t)

def wavePhi_init(x,t):
    return wavePhi(x,t)

def waveVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t))

def waveVF_init(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi_init(x,t))

def twpflowVelocity_u(x,t):
#    waterspeed = waveVelocity_u(x,t)
#    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
#    return H*windspeed_u + (1.0-H)*waterspeed
    return Um

def twpflowVelocity_v(x,t):
#    waterspeed = waveVelocity_v(x,t)
#    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
#    return H*windspeed_v+(1.0-H)*waterspeed
    return 0.0

def twpflowVelocity_w(x,t):
#    waterspeed = waveVelocity_w(x,t)
#    H = smoothedHeaviside(epsFact_consrv_heaviside*he,wavePhi(x,t)-epsFact_consrv_heaviside*he)
#    return H*windspeed_w+(1.0-H)*waterspeed
    return 0.0

def twpflowVelocity_u_init(x,t):
    return twpflowVelocity_u(x,t)

def twpflowVelocity_v_init(x,t):
    return twpflowVelocity_v(x,t)

def twpflowVelocity_w_init(x,t):
    return twpflowVelocity_w(x,t)

def twpflowFlux(x,t):
    return -twpflowVelocity_u(x,t)

def outflowVF(x,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he,x[2] - outflowHeight)

def twpflowPressure(x,t):
    #return min(L[2] - x[2],L[2]-waterLevel)*rho_1*(-g[2]) + max(waterLevel - x[2],0.0)*rho_0*(-g[2])
    p_L = L[2]*rho_1*g[2]
    phi_L = wavePhi((x[0],x[1],L[2]),t)
    phi = wavePhi(x,t)
    return p_L - g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                          -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))
def twpflowPressure_init(x,t):
    #return min(L[2] - x[2],L[2]-waterLevel)*rho_1*(-g[2]) + max(waterLevel - x[2],0.0)*rho_0*(-g[2])
    p_L = L[2]*rho_1*g[2]
    phi_L = L[2] - inflowHeightMean
    phi = x[2] - inflowHeightMean
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

def outflowPressure(x,t):
    #return min(L[2] - x[2],L[2]-waterLevel)*rho_1*(-g[2]) + max(waterLevel - x[2],0.0)*rho_0*(-g[2])
    p_L = L[2]*rho_1*g[2]
    phi_L = L[2] - outflowHeight
    phi = x[2] - outflowHeight
    return p_L -g[2]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))
