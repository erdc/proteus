cimport cython
import cython

cdef class BC_Base:
    cdef public:
        object Shape
        str name
        str BC_type
        double[:] _b_or
        object ct
    # cpdef void newGlobalBC(BC_Base cls, str name, object default_value)
    # cpdef void getContext(BC_Base cls, object context=*)

ctypedef double (*cpp_uOfXT) (BoundaryCondition, double[:], double)
cdef class BoundaryCondition:
    cdef cpp_uOfXT uuOfXT
    cdef public:
        object uOfXT
    cpdef void resetBC(self)
    # cpdef void setConstantBC(self, double value)
    # cpdef void setLinearBC(self, double a0, double a1, int i)
