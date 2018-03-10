from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent
from proteus.mprans import SpatialTools as st
from proteus import Context

ct = Context.Options([
    ("T", 4.0, "Time interval [0, T]"),
    ("Refinement",16, "refinement"),
    ("onlySaveFinalSolution",False,"Only save the final solution"),
    ("vspaceOrder",1,"FE space for velocity"),
    ("pspaceOrder",1,"FE space for pressure"),
    ("parallel",False,"Use parallel or not"),
    ("use_sbm",False,"use sbm instead of imb"),
    ("is_structured",False,"use structured mesh or not"),
    ("free_x",(0,1,0),"which direction is allowed"),
    ("free_r",(0,0,1),"which direction is allowed"),
    ("solid_dt",0.00001,"time step for solid motion")
], mutable=True)


if ct.use_sbm:
    USE_SBM=1
else:
    USE_SBM=0

#  Discretization -- input options
#Refinement = 20#45min on a single core for spaceOrder=1, useHex=False
Refinement = ct.Refinement
sedimentDynamics=False
genMesh = True
movingDomain = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = True


parallel = ct.parallel
if parallel:
    usePETSc = True
    useSuperlu=False
else:
    usePETSc = False
parallelPartitioningType = MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1

timeDiscretization = 'vbdf'#vbdf'#'vbdf'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = ct.vspaceOrder
pspaceOrder = ct.pspaceOrder
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 0.0
useOnlyVF = False
useRANS = 0  # 0 -- None
             # 1 -- K-Epsilon
             # 2 -- K-Omega
openTop=True
fl_H = 0.41
# Input checks
if spaceOrder not in [1, 2]:
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
    hFactor = 1.0
    if useHex:
        basis = C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 2)
    else:
        basis = C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 5)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 5)

if pspaceOrder == 1:
    if useHex:
        pbasis = C0_AffineLinearOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineLinearOnSimplexWithNodalBasis
elif pspaceOrder == 2:
    if useHex:
        pbasis = C0_AffineLagrangeOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineQuadraticOnSimplexWithNodalBasis

# Domain and mesh
#L = (0.584,0.350)
x0 = (0, 0)
L  = (2, 5)
he = 0.04
#he*=0.5
#he*=0.5
#he*=0.5
weak_bc_penalty_constant = 100.0
nLevels = 1
# parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured = ct.is_structured



# new style PointGauges
# pointGauges = PointGauges(gauges = ((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))), (('p',), ((0.5, 0.5, 0),))),
#                           activeTime=(0, 0.5),
#                           sampleRate=0,
#                           fileName='combined_gauge_0_0.5_sample_all.csv')

# lineGauges = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825': ((0.495, 0.0, 0.0), (0.495, 1.8, 0.0))}, linePoints=20)
# #'lineGauge_x/H=1.653':((0.99,0.0,0.0),(0.99,1.8,0.0))
# lineGauges_phi = LineGauges_phi(lineGauges.endpoints, linePoints=20)

if useHex:
    nnx = 4 * Refinement + 1
    nny = 2 * Refinement + 1
    hex = True
    domain = Domain.RectangularDomain(L)
else:
    boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    if structured:
        nnx = 2 * Refinement + 1
        nny = 5 * Refinement + 1
    else:
        vertices = [[x0[0], x0[1]],  #0
                    [x0[0]+L[0], x0[1]],  #1
                    [x0[0]+L[0], x0[1]+L[1]],  #2
                    [x0[0], x0[1]+L[1]],  #3
                    ]
        vertexFlags = [boundaryTags['bottom'],
                       boundaryTags['bottom'],
                       boundaryTags['top'],
                       boundaryTags['top'],
                       # the interior vertices should be flaged to 0
                       ]

        segments = [[0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0],
                    #Interior segments
                    ]
        segmentFlags = [boundaryTags['bottom'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['left']]

        regions = [[x0[0]+0.5*L[0],x0[1]+0.5*L[1]]]
        regionFlags = [1]
        regionConstraints=[0.5*he**2]
        #        for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
        #            vertices.append(gaugeCoordinates)
        #            vertexFlags.append(pointGauges.flags[gaugeName])

        # for gaugeName, gaugeLines in lineGauges.linepoints.iteritems():
        #     for gaugeCoordinates in gaugeLines:
        #         vertices.append(gaugeCoordinates)
        #         vertexFlags.append(lineGauges.flags[gaugeName])
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags,
                                                      regionConstraints=regionConstraints)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        #triangleOptions = "VApq30Dena%8.8f" % ((he ** 2) / 2.0,)
        triangleOptions = "VApq30Dena"

        logEvent("""Mesh generated using: tetgen -%s %s""" % (triangleOptions, domain.polyfile + ".poly"))
# Time stepping
T=ct.T
dt_fixed = 0.005#0.03
dt_init = 0.0025#min(0.1*dt_fixed,0.001)
runCFL=0.33
nDTout = int(round(T/dt_fixed))
if dt_init<dt_fixed:
    tnList = [0.0,dt_init]+[i*dt_fixed for i in range(1,nDTout+1)]
else:
    tnList = [i*dt_fixed for i in range(0,nDTout+1)]

if ct.onlySaveFinalSolution == True:
    tnList = [0.0,dt_init,ct.T]
# Numerical parameters
ns_forceStrongDirichlet = True
ns_sed_forceStrongDirichlet = False
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
    kappa_shockCapturingFactor = 0.25
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.25
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ns_sed_shockCapturingFactor = 0.9
    ns_sed_lag_shockCapturing = True
    ns_sed_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    vos_shockCapturingFactor = 0.9
    vos_lag_shockCapturing = True
    vos_sc_uref = 1.0
    vos_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_vos = epsFact_consrv_heaviside = epsFact_consrv_dirac = \
        epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 0.1
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0

ns_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
ns_sed_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
vos_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.05 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
phi_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
pressure_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)

#turbulence
ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
# Water
rho_0 = 1.0
nu_0 = 0.01

# Air
rho_1 = rho_0#1.205
nu_1 = nu_0#1.500e-5

# Sediment
rho_s = rho_0
nu_s = 10000.0*nu_0
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0, -981]

# Initial condition
waterLine_x = 0.75
waterLine_z = 1.6

def quarter_circle(center, radius, p_nb, angle, angle0=0., v_start=0.):
    vertices = []
    segments = []
    for i in range(p_nb+1):
        x = radius*np.sin(angle0+angle*float(i)/(p_nb))
        y = radius*np.cos(angle0+angle*float(i)/(p_nb))
        vertices += [[center[0]+x, center[1]+y]]
        if i > 0:
            segments += [[v_start+(i-1), v_start+i]]
        elif i == p_nb-1:
            segments += [[v_start+i, v_start]]
    return vertices, segments

def get_a_cylinder_shape(_d,_c, _r, _bdTag, _b=(0.0,0.0,0.0), _a=0.0):#create a cylinder with center _c and radius _r, retate about _b angle _a, with boundary flag _f
    p_nb = int((np.pi*2*_r)/(he))
    v, s = quarter_circle(center=[0.,0.], radius=_r, p_nb=p_nb,angle=2*np.pi, angle0=0., v_start=0.)
    vertices = []
    vertexFlags = []
    segments = []
    segmentFlags = []
    center = [0., 0.]
    flag = 1
    v_start = 0
    vertices += v[:-1]
    vertexFlags += [1]*len(vertices)
    segments += s[:-1]+[[len(vertices)-1, 0]]
    segmentFlags += [1]*len(segments)
    segments[-1][1] = 0  # last segment links to vertex 0
    boundaryTags = {_bdTag: 1}
    caisson = st.CustomShape(_d, barycenter=_b,
                            vertices=vertices, vertexFlags=vertexFlags,
                            segments=segments, segmentFlags=segmentFlags,
                            boundaryTags=boundaryTags)
    facet = []
    for i, vert in enumerate(caisson.vertices):
        facet += [i]
    caisson.facets = np.array([[facet]])
    caisson.facetFlags = np.array([1])
    caisson.regionFlags = np.array([1])
    
    ang = _a
    caisson.rotation_init = np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
    caisson.rotate(ang, pivot=caisson.barycenter)
    
    caisson.setHoles([[0., 0.]])
    caisson.holes_ind = np.array([0])
    caisson.translate([_c[0], _c[1]])
    
    for bc in caisson.BC_list:
        bc.setNoSlip()
    return caisson

def add_shape(_shape, _g, _M, _Iz, _filename):
    from proteus import AuxiliaryVariables
    class BodyAddedMass(AuxiliaryVariables.AV_base):
        def __init__(self):
            AuxiliaryVariables.AV_base.__init__(self)

        def attachModel(self, model, ar):
            self.model=model
            return self

        def calculate(self):
            pass
    from proteus.mprans import BodyDynamics as bd
    body = bd.RigidBody(shape=_shape,g=_g)
#     system = body  # so it is added in twp_navier_stokes_p.py
    bodyAddedMass = BodyAddedMass()
    body.bodyAddedMass = body.ProtChAddedMass = bodyAddedMass
    body.setMass(_M)
    free_x = ct.free_x
    free_r = ct.free_r
    body.setConstraints(free_x=free_x,
                        free_r=free_r)
    body.It = np.array([[1., 0., 0.],
                        [0., 1., 0.],
                        [0., 0., _Iz]])
    filename=_filename
    body.setRecordValues(filename=filename, all_values=True)
    body.coords_system = _shape.coords_system  # hack
    body.last_Aij = np.zeros((6,6),'d')
    body.last_a = np.zeros((6,),'d')
    body.a = np.zeros((6,),'d')
    body.last_Omega = np.zeros((3,3),'d')
    body.last_velocity = np.zeros((3,),'d')
    body.last_Q = np.eye(3)
    body.last_h = np.zeros((3,),'d')
    body.last_mom = np.zeros((6,),'d')
    body.last_u = np.zeros((18,),'d')
    body.last_t = 0.0
    body.t = 0.0
    body.h = np.zeros((3,),'d')
    body.velocity = np.zeros((3,),'d')
    body.mom = np.zeros((6,),'d')
    body.Omega = np.zeros((3,3),'d')
    body.Q = np.eye(3,3)
    body.u = np.zeros((18,),'d')
    body.u[9:18]= body.Q.flatten()
    body.last_u[:] = body.u
    body.last_Q[:] = body.Q
    body.FT = np.zeros((6,),'d')
    body.last_FT = np.zeros((6,),'d')
    body.free_dof = np.zeros((6,),'d')
    body.free_dof[:3] = free_x
    body.free_dof[3:] = free_r
    # cek -- I  did this with a different hack, see dummy AddedMassBody above
    # body.Aij = np.zeros((6,6))

    # # we create a ProtChAddedMass auxiliary variable to add to added_mass_n.py
    # # the argument here should be a ProtChSystem, so we give it an empty system (hack)
    # # this is just to get access to the AddedMass model from the body functions
    # system = crb.ProtChSystem(g)
    # body.ProtChAddedMass = crb.ProtChAddedMass(system)

    def step(dt):
        print "yy1:", body.F
        print "yy2:", body.M
        print "yy3:", body.position
        
        from math import ceil
        logEvent("Barycenter "+str(body.barycenter))
        n = max(1.0,ceil(dt/ct.solid_dt))
        DT=dt/n
        def F(u,theta):
            """The residual and Jacobian for discrete 6DOF motion with added mass"""
            v = u[:3]
            omega = u[3:6]
            h = u[6:9]
            Q = u[9:18].reshape(3,3)
            Omega = np.array([[      0.0, -omega[2],  omega[1]],
                              [ omega[2],       0.0, -omega[0]],
                              [-omega[1],  omega[0],      0.0]])
            I = np.matmul(np.matmul(Q, body.It), Q.transpose())
            body.Aij = np.zeros((6,6),'d')
#             if opts.addedMass:
#                 for i in range(1,5):#number of rigid body facets
#                     body.Aij += body.bodyAddedMass.model.levelModelList[-1].Aij[i]
            avg_Aij=False
            if avg_Aij:
                M = body.Aij*theta + body.last_Aij*(1-theta)
            else:
                M = body.Aij.copy()
            for i in range(6):
                for j in range(6):
                    M[i,j]*=body.free_dof[j]#only allow j accelerations to contribute to i force balance if j is free
                    M[j,i]*=body.free_dof[j]#only allow j added mass forces if j is free
                    body.Aij[i,j]*=body.free_dof[j]#only allow j accelerations to contribute to i force balance if j is free
                    body.Aij[j,i]*=body.free_dof[j]#only allow j added mass forces if j is free
            body.FT[:3] = body.model.levelModelList[-1].coefficients.particle_netForces[0]#body.F 
            body.FT[3:] = body.model.levelModelList[-1].coefficients.particle_netMoments[0]#body.M
            body.FT[1] += body.mass * _g[1]
            
            print "yy4:", body.FT[:3]
            print "yy5:", body.FT[3:]
            for i in range(3):
                M[i, i] += body.mass
                for j in range(3):
                    M[3+i, 3+j] += I[i, j]
            r = np.zeros((18,),'d')
            BE=True
            CF=1
            if BE:
                r[:6] = np.matmul(M, u[:6]) - np.matmul(body.Aij, body.last_u[:6]) - body.last_mom - DT*body.FT - CF*np.matmul(body.Aij,body.last_a)
                r[6:9] = h - body.last_h - DT*v
                rQ = Q - body.last_Q - DT*np.matmul(Omega,Q)
            else:
                r[:6] = np.matmul(M, u[:6]) - np.matmul(body.Aij, body.last_u[:6]) - body.last_mom - DT*(body.FT*theta+body.last_FT*(1.0-theta)) - CF*np.matmul(body.Aij,body.last_a)
                r[6:9] = h - body.last_h - DT*0.5*(v + body.last_velocity)
                rQ = Q - body.last_Q - DT*0.25*np.matmul((Omega + body.last_Omega),(Q+body.last_Q))
            r[9:18] = rQ.flatten()
            J = np.zeros((18,18),'d')
            J[:6,:6] = M
            #neglecting 0:6 dependence on Q
            for i in range(3):
                if BE:
                    J[6+i,i] = -DT
                else:
                    J[6+i,i] = -DT*0.5
                J[6+i,6+i] = 1.0
            for i in range(9):
                J[9+i,9+i] = 1.0
            for i in range(3):
                for j in range(3):
                    if BE:
                        J[9+i*3+j, 9+i+j*3] -= DT*Omega[i,j]
                    else:
                        J[9+i*3+j, 9+i+j*3] -= DT*0.25*(Omega+body.last_Omega)[i,j]
            #neglecting 9:18 dependence on omega
            body.Omega[:] = Omega
            body.velocity[:] = v
            body.Q[:] = Q
            body.h[:] = h
            body.u[:] = u
            body.mom[:3] = body.mass*u[:3]
            body.mom[3:6] = np.matmul(I,u[3:6])
            body.a[:] = u[:6] - body.last_u[:6]
            return r, J
        nd = body.nd
        Q_start=body.Q.copy()
        h_start=body.last_h.copy()
        for i in range(int(n)):
            theta = (i+1)*DT/dt
            body.t = body.last_t + DT
            logEvent("6DOF theta "+`theta`)
            u = np.zeros((18,),'d')
            u[:] = body.last_u
            r = np.zeros((18,),'d')
            r,J = F(u,theta)
            its=0
            maxits=100
            while ((its==0 or np.linalg.norm(r) > 1.0e-8) and its < maxits):
                u -= np.linalg.solve(J,r)
                r,J = F(u,theta)
                its+=1
            logEvent("6DOF its "+`its`)
            logEvent("6DOF res "+`np.linalg.norm(r)`)
            if its==maxits:
                import pdb
                pdb.set_trace()
                print its,np.linalg.norm(r)
            body.last_Aij[:]=body.Aij
            body.last_FT[:] = body.FT
            body.last_Omega[:] = body.Omega
            body.last_velocity[:] = body.velocity
            body.last_Q[:] = body.Q
            body.last_h[:] = body.h
            body.last_u[:] = body.u
            body.last_mom[:] = body.mom
            body.last_t = body.t
        body.last_a[:] = body.a
        # translate and rotate
        body.h -= h_start
        body.last_h[:] = 0.0
        body.last_position[:] = body.position
        #body.rotation_matrix[:] = np.linalg.solve(Q_start,body.Q)
#         body.rotation_matrix[:] = np.matmul(np.linalg.inv(body.Q),Q_start)
#         body.rotation_euler[2] -= math.asin(body.rotation_matrix[1,0])#cek hack
        body.Shape.translate(body.h[:nd])
        body.barycenter[:] = body.Shape.barycenter
        body.position[:] = body.Shape.barycenter
        body.model.levelModelList[-1].coefficients.particle_centroids[0][:] = body.Shape.barycenter
        print "yy6:", body.position
        #logEvent("6DOF its = " + `its` + " residual = "+`r`)
        logEvent("6DOF time "+`body.t`)
        logEvent("6DOF DT "+`DT`)
        logEvent("6DOF n "+`n`)
        logEvent("6DOF force "+`body.FT[1]`)
        logEvent("displacement, h = "+`body.h`)
        logEvent("rotation, Q = "+`body.Q`)
        logEvent("velocity, v = "+`body.velocity`)
        logEvent("angular acceleration matrix, Omega = "+`body.Omega`)

    body.step = step
    #body.scheme = 'Forward_Euler'
    return body


R = 0.125
Rho_S=1.5
M = np.pi*R*R*Rho_S
Iz= M*0.5*R*R
C1 = (1.0,4.0)
C2 = (1.0,3.5)

caisson1 = get_a_cylinder_shape(domain,C1, R, 'c1')
cylinder1 = add_shape(caisson1,np.array([0,-981,0]),M,Iz,'c1')

def particle_sdf_1(t, x):
    cx = cylinder1.barycenter[0]
    cy = cylinder1.barycenter[1]
    r = math.sqrt( (x[0]-cx)**2 + (x[1]-cy)**2)
    n = ((x[0]-cx)/(r+1e-10),(x[1]-cy)/(r+1e-10))
    
    return  r-R,n

def particle_vel_1(t, x):
    r = np.array([x[0]-cylinder1.barycenter[0], x[1]-cylinder1.barycenter[1],0])#use barycenter not position since the first step error
    u = np.dot(cylinder1.Omega, r) + cylinder1.velocity
    #return cylinder1.velocity[:2]
    return u[:2]
