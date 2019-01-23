cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)
from cython.operator cimport dereference as deref
cimport ChronoHeaders as ch

cdef class ChBodyAddedMass:
    cdef shared_ptr[ch.ChBodyAddedMass] sharedptr
    cdef shared_ptr[ch.ChBody] sharedptr_chbody
    cdef ch.ChBodyAddedMass * thisptr
    cdef ch.ChBody * bodyptr
    cdef void SetMfullmass(self, ch.ChMatrixDynamic Mfullmass_in)
    cdef void SetInvMfullmass(self, ch.ChMatrixDynamic Mfullmass_in)
    cdef public:
        object ChBodySWIG
