# distutils: language = c++
"""
This file defines python object using a syntax as close as possible to the
Chrono syntax. This creates and gives access to C++ objects.
"""

cimport numpy as np
import numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)
from cython.operator cimport dereference as deref
cimport ChronoHeaders as ch
import ChronoEngine_python_core as chrono
from libc.stdlib cimport malloc, free


cdef extern from "swigpyobject.h":
    ctypedef struct SwigPyObject:
        void *ptr


cdef class ChBodyAddedMass:
    """Cython class for ChBodyAddedMass
    (!) Uses shared_ptr
    """

    def __cinit__(self):
        # make shared_ptr object (C++ syntax for Chrono)
        self.sharedptr = make_shared[ch.ChBodyAddedMass]()
        self.sharedptr_chbody = <shared_ptr[ch.ChBody]> self.sharedptr
        # get a raw ptr for SWIG
        self.thisptr = self.sharedptr.get()
        self.bodyptr = self.sharedptr_chbody.get()
        # create SWIG ChBody
        self.ChBodySWIG = chrono.ChBody()
        # delete? object pointed to by SWIG
        self.ChBodySWIG.this.disown()
        # point to new object (base of ChBodyAddedMass: ChBody)
        cdef SwigPyObject *swig_obj = <SwigPyObject*>self.ChBodySWIG.this
        swig_obj.ptr = <ch.ChBody*?> &self.bodyptr

    cdef void SetMfullmass(self, ch.ChMatrixDynamic Mfullmass_in):
        self.thisptr.SetMfullmass(Mfullmass_in)

    cdef void SetInvMfullmass(self, ch.ChMatrixDynamic inv_Mfullmass_in):
        self.thisptr.SetInvMfullmass(inv_Mfullmass_in)

    cpdef void SetMass(self, double newmass):
        self.thisptr.SetMass(newmass)

    cpdef double GetMass(self):
        return self.thisptr.GetMass()

    cpdef double GetPos(self):
        return self.thisptr.GetPos().x()
