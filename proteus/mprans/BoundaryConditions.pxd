cimport cython
cimport numpy as np
from proteus.BoundaryConditions cimport (BC_Base,
                                         BoundaryCondition)
ctypedef double[:,:,:] double_memview3
ctypedef double[:,:] double_memview2
ctypedef double[:] double_memview1
ctypedef int[:] int_memview1

# ctypedef double (*cpp_uOfXT) (BC_RANS, double[:], double)

cdef class BC_RANS(BC_Base):
    cdef double[:] zero_array
    cdef double f
    cdef __cppClass_WavesCharacteristics waves
    cdef object body
    cdef double[:] wind_speed
    cdef double[:, :] body_python_rot_matrix
    cdef double[:] body_python_last_pos
    cdef double[:] body_python_h
    cdef public:
        # dirichlet
        cdef BoundaryCondition p_dirichlet
        cdef BoundaryCondition u_dirichlet
        cdef BoundaryCondition v_dirichlet
        cdef BoundaryCondition w_dirichlet
        cdef BoundaryCondition vof_dirichlet
        cdef BoundaryCondition k_dirichlet
        cdef BoundaryCondition dissipation_dirichlet
        # advective
        cdef BoundaryCondition p_advective
        cdef BoundaryCondition u_advective
        cdef BoundaryCondition v_advective
        cdef BoundaryCondition w_advective
        cdef BoundaryCondition vof_advective
        cdef BoundaryCondition k_advective
        cdef BoundaryCondition dissipation_advective
        #diffusive
        cdef BoundaryCondition u_diffusive
        cdef BoundaryCondition v_diffusive
        cdef BoundaryCondition w_diffusive
        cdef BoundaryCondition k_diffusive
        cdef BoundaryCondition dissipation_diffusive
        # mesh
        cdef BoundaryCondition hx_dirichlet
        cdef BoundaryCondition hy_dirichlet
        cdef BoundaryCondition hz_dirichlet
        cdef BoundaryCondition u_stress
        cdef BoundaryCondition v_stress
        cdef BoundaryCondition w_stress
        # functions
        # cpdef void reset(self)
        # cpdef void setNonMaterial(self)
        # cpdef void setTank(self)
        # cpdef void setFixedNodes(self)
        # cpdef void setNoSlip(self)
        # cpdef void setFreeSlip(self)
        # cpdef void setAtmosphere(self, double[:] orientation=*, double vof_air=*)
        # cpdef void setMoveMesh(self, double[:] last_pos, double[:] h=*, double[:,:] rot_matrix=*)
        cpdef double* __cpp_MoveMesh_h(self, double[:] x, double t)
        cpdef double __cpp_MoveMesh_hx(self, double[:] x, double t)
        cpdef double __cpp_MoveMesh_hy(self, double[:] x, double t)
        cpdef double __cpp_MoveMesh_hz(self, double[:] x, double t)
#         @cython.locals(waveHeight=double, wavePhi=double,
#                        water_speed=double_memview1, x_max=double_memview1,
#                        u=cython.p_double, H=double)
        cdef double __cpp_UnsteadyTwoPhaseVelocityInlet_u_dirichlet(self, double[:] x, double t)
        cdef double __cpp_UnsteadyTwoPhaseVelocityInlet_v_dirichlet(self, double[:] x, double t)
        cdef double __cpp_UnsteadyTwoPhaseVelocityInlet_w_dirichlet(self, double[:] x, double t)
        cdef double __cpp_UnsteadyTwoPhaseVelocityInlet_vof_dirichlet(self, double[:] x, double t)
        cdef double __cpp_UnsteadyTwoPhaseVelocityInlet_p_advective(self, double[:] x, double t)

ctypedef double[:] (*cfvel) (double[3], double)  # pointer to velocity function
ctypedef double (*cfeta) (double[3], double)  # pointer to eta function
ctypedef double[:] (*cfvelrel) (RelaxationZone, double[3], double)  # pointer to velocity function of RelaxationZone class
ctypedef double (*cfphirel) (RelaxationZone, double[3])  # pointer to phi function of RelaxationZone class
ctypedef np.float64_t float64_t
ctypedef np.int64_t int64_t

cdef class RelaxationZone:
    cdef int vert_axis
    cdef double[:] wind_speed
    cdef double[:] u_calc  # calculated velocity
    # wave characteristics (if any)
    cdef cfvelrel uu
    cdef cfphirel phi
    cdef double mwl
    cdef cfeta eta
    cdef double waveHeight
    cdef double wavePhi
    cdef double waterSpeed
    # for smoothing
    cdef double he
    cdef double ecH
    cdef double H
    cdef int nd
    cdef double[:] zero_vel
    cdef __cppClass_WavesCharacteristics waves
    cdef public:
        object Shape
        str zone_type
        double dragAlpha
        double dragBeta
        double porosity
        double epsFact_solid
        double[:] center
        double[:] orientation
    cpdef void calculate_init(self)
    cdef double[:] calculate_vel(self, cython.p_double x, double t)
    @cython.locals(xx=double_memview1, waveHeight=double, wavePhi=double,
                   waterSpeed=double_memview1, x_max=double_memview1, H=double,
                   u=cython.double[3], o1=double, o2=double, o3=double)
    cdef double[:] __cpp_calculate_vel_wave(self, double* x, double t)
    cdef double[:] __cpp_calculate_vel_zero(self, double* x, double t)
    cdef double calculate_phi(self, double[3] x)
    cdef double __cpp_calculate_phi(self, double[3] x)
    @cython.locals(d1=double, d2=double, d3=double, phi=double,
                   o1=double, o2=double, o3=double)
    cdef double __cpp_calculate_phi_porous(self, double[3] x)

cdef class RelaxationZoneWaveGenerator:
    cdef int nd  # dimension
    cdef int max_flag  # maximum region flag of relaxation zones (initialised in calculate_init)
    cdef RelaxationZone[:] zones_array  # zones array for fast access
    cdef public:
        dict zones  # zones dictionary
        object model  # model attached to zone
        object ar  #
    cpdef RelaxationZoneWaveGenerator attachModel(self, dict, int)  # needed
    cpdef void attachAuxiliaryVariables(self, dict)  # needed
    @cython.locals(zones=np.ndarray)
    cpdef void calculate_init(self)
    cpdef void calculate(self)
    @cython.locals(m=object, zone=RelaxationZone, mType=int, nE=int, nk=int,
                    qx=double_memview3, x=cython.double[3], t=double, nl=int,
                    q_phi_solid=double_memview2,  phi=double,
                   q_velocity_solid=double_memview3, u=double_memview1, mTypes=int_memview1)
    cdef void __cpp_iterate(self)  # main iteration loop


cdef class __cppClass_WavesCharacteristics:
    cdef int vert_axis  # vertical axis
    # relaxation zone info (if any)
    cdef double[:] center  # center of zone
    cdef double[:] orientation  # orientation of zone
    # zero array
    cdef double[:] zero_vel
    cdef double[:] _b_or
    cdef double[:] wind_speed
    cdef public:
        # wave class from WaveTools
        object WT
    @cython.locals(waveHeight=double, wavePhi=double,
                   waterSpeed=double_memview1, x_max=double_memview1, H=double,
                   o1=double, o2=double, o3=double)
    cdef double[:]  __cpp_calculate_velocity(self, double[3] x, double t)
    @cython.locals(_b_or=double_memview1, ux=double_memview1)
    cdef double __cpp_calculate_pressure(self, double[3] x, double t)
    @cython.locals(level=double)
    cdef double __cpp_calculate_vof(self, double[3] x, double t)
    cdef double __cpp_calculate_phi(self, double[3] x)

cdef double* __x_to_cpp(double[:])
