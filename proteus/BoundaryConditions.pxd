cimport cython
import cython

cdef class BC_Base:
    cdef public:
        object Shape
        str name
        str BC_type
        double[:] _b_or
        int _b_i
    cpdef void newGlobalBC(BC_Base cls, str name, object default_value)
    cpdef void getContext(BC_Base cls, object context=*)

ctypedef double (*BC_uOfXT) (double[3], double)
cdef class BoundaryCondition:
    cdef BC_uOfXT uOfXT


@cython.locals(testt=double, a=double, b=double)
cdef double test0(double[:], double)
