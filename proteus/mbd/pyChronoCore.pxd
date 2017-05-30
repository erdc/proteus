cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)
from cython.operator cimport dereference as deref
cimport ChronoHeaders as ch

cdef class ChVector:
    cdef ch.ChVector cppobj
    cpdef double x(self)
    cpdef double y(self)
    cpdef double z(self)

cdef class ChQuaternion:
    cdef ch.ChQuaternion cppobj
    cpdef double e0(self)
    cpdef double e1(self)
    cpdef double e2(self)
    cpdef double e3(self)

cdef class ChFrame:
    cdef shared_ptr[ch.ChFrame] sharedptr_chframe
    cpdef np.ndarray GetPos(self)
    cpdef np.ndarray GetRot(self)
    cpdef void SetRot(self, ChQuaternion rot)
    cpdef void SetPos(self, ChVector mpos)
    cpdef np.ndarray GetA(self)
    cpdef np.ndarray GetRotAxis(self)
    cpdef double GetRotAngle(self)

cdef class ChFrameMoving(ChFrame):
    cdef shared_ptr[ch.ChFrameMoving] sharedptr_chframemoving
    cpdef np.ndarray GetPos_dt(self)
    cpdef np.ndarray GetPos_dtdt(self)
    cpdef np.ndarray GetRot_dt(self)
    cpdef np.ndarray GetRot_dtdt(self)
    cpdef np.ndarray GetWvel_loc(self)
    cpdef np.ndarray GetWacc_loc(self)

cdef class ChBodyFrame(ChFrameMoving):
    cdef shared_ptr[ch.ChBodyFrame] sharedptr_chbodyframe

cdef class ChBody(ChBodyFrame):
    cdef shared_ptr[ch.ChBody] sharedptr_chbody
    cpdef void SetBodyFixed(self, bool state)
    cpdef void SetMaterialSurface(self, ChMaterialSurfaceSMC mat)
    cpdef void SetInertiaXX(self, ChVector iner)
    cpdef void SetInertiaXY(self, ChVector iner)
    cpdef void SetMass(self, double newmass)
    cpdef double GetMass(self)

cdef class ChBodyEasyBox(ChBody):
    cdef shared_ptr[ch.ChBodyEasyBox] sharedptr_chbodyeasybox

cdef class ChMaterialSurfaceSMC:
    cdef shared_ptr[ch.ChMaterialSurfaceSMC] sharedptr
    cpdef void SetYoungModulus(self, float val)
    cpdef void SetPoissonRatio(self, float val)
    cpdef void SetSfriction(self, float val)
    cpdef void SetKfriction(self, float val)
    cpdef void SetFriction(self, float val)
    cpdef void SetRestitution(self, float val)
    cpdef void SetAdhesion(self, float val)
    cpdef void SetAdhesionMultDMT(self, float val)
    cpdef void SetKn(self, float val)
    cpdef void SetKt(self, float val)
    cpdef void SetGn(self, float val)
    cpdef void SetGt(self, float val)


cdef class ChContactSurfaceNodeCloud:
    cdef shared_ptr[ch.ChContactSurfaceNodeCloud] sharedptr
    cpdef void AddAllNodes(self, double point_radius)



cdef np.ndarray ChVector_to_npArray(ch.ChVector &myvec)
cdef np.ndarray ChQuaternion_to_npArray(ch.ChQuaternion &quat)
cdef np.ndarray ChMatrix33_to_npArray(ch.ChMatrix33 &mat)
cdef ch.ChVector npArray_to_ChVector(np.ndarray arr)
cdef ch.ChQuaternion npArray_to_ChQuaternion(np.ndarray arr)

# cdef class ChLinkPointPoint:
#     cdef shared_ptr[ch.ChLinkPointPoint] sharedptr
#     def __cinit__(self):
#         self.sharedptr = make_shared[ch.ChLinkPointPoint]()
#     def Initialize(shared_ptr[ChNodeFEAxyz] anodeA, shared_ptr[ChNodeFEAxyz])
#         deref(self.sharedptr).Initialize(<shared_ptr[ch.ChNodeFEAxyz]> anodeA, <shared_ptr[ch.ChNodeFEAxyz]> anodeB)
