cimport cython
cimport numpy as np
from proteus.BoundaryConditions cimport (BC_Base,
                                         BoundaryCondition)
ctypedef double[:,:,:] double_memview3
ctypedef double[:,:] double_memview2
ctypedef double[:] double_memview1
ctypedef int[:] int_memview1


cdef class BC_RANS(BC_Base):
    cdef double[:] zero_array
    cdef double f
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
        cdef double u_stress
        cdef double v_stress
        cdef double w_stress
        # functions
        cpdef void reset(self)
        cpdef void setNonMaterial(self)
        cpdef void setTank(self)
        cpdef void setFixedNodes(self)
        cpdef void setNoSlip(self)
        cpdef void setFreeSlip(self)
        cpdef void setAtmosphere(self, double[:] orientation=*, double vof_air=*)
        cpdef void setMoveMesh(self, double[:] last_pos, double[:] h=*, double[:,:] rot_matrix=*)

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
    cdef cfvel u
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
    cdef public:
        object Shape
        str zone_type
        double dragAlpha
        double dragBeta
        double porosity
        double epsFact_solid
        object waves
        double[:] center
        double[:] orientation
    cpdef void calculate_init(self)
    cdef double[:] calculate_vel(self, cython.p_double x, double t)
    @cython.locals(xx=double_memview1, waveHeight=double, wavePhi=double,
                   waterSpeed=double_memview1, x_max=double_memview1, H=double, u=double[3],
                   o1=double, o2=double, o3=double)
    cdef double[:] __cpp_calculate_vel_wave(self, cython.p_double x, double t)
    cdef double[:] __cpp_calculate_vel_zero(self, cython.p_double x, double t)
    cdef double calculate_phi(self, cython.p_double x)
    @cython.locals(d1=double, d2=double, d3=double, phi=double,
                   o1=double, o2=double, o3=double)
    cdef double __cpp_calculate_phi(self, cython.p_double x)
    cdef double __cpp_calculate_phi_porous(self, cython.p_double x)

cdef class RelaxationZoneWaveGenerator:
    cdef int nd
    cdef public:
        dict zones
        object model
        object ar
    cpdef RelaxationZoneWaveGenerator attachModel(self, dict, int)
    cpdef void attachAuxiliaryVariables(self, dict)
    @cython.locals(zones=np.ndarray)
    cpdef void calculate_init(self)
    cpdef void calculate(self)
    @cython.locals(m=object, zone=RelaxationZone, mType=int, nE=int, nk=int,
                    qx=double_memview3, x=cython.p_double, t=double, nl=int,
                    q_phi_solid=double_memview2,  phi=double,
                   q_velocity_solid=double_memview3, u=double_memview1, mTypes=int_memview1)
    cdef void __cpp_iterate(self)

