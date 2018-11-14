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


cdef np.ndarray ChVector_to_npArray(ch.ChVector &myvec):
    return np.array([myvec.x(), myvec.y(), myvec.z()])

cdef np.ndarray ChQuaternion_to_npArray(ch.ChQuaternion &quat):
    return np.array([quat.e0(), quat.e1(), quat.e2(), quat.e3()])

cdef np.ndarray ChMatrix33_to_npArray(ch.ChMatrix33 &mat):
    return np.array([[mat.Get_A_Xaxis().x(), mat.Get_A_Xaxis().y(), mat.Get_A_Xaxis().z()],
                     [mat.Get_A_Yaxis().x(), mat.Get_A_Yaxis().y(), mat.Get_A_Yaxis().z()],
                     [mat.Get_A_Zaxis().x(), mat.Get_A_Zaxis().y(), mat.Get_A_Zaxis().z()]])

cdef np.ndarray ConstChMatrix33_to_npArray(const ch.ChMatrix33 &mat):
    return np.array([[mat.Get_A_Xaxis().x(), mat.Get_A_Xaxis().y(), mat.Get_A_Xaxis().z()],
                     [mat.Get_A_Yaxis().x(), mat.Get_A_Yaxis().y(), mat.Get_A_Yaxis().z()],
                     [mat.Get_A_Zaxis().x(), mat.Get_A_Zaxis().y(), mat.Get_A_Zaxis().z()]])

cdef ch.ChVector npArray_to_ChVector(np.ndarray arr):
    cdef ch.ChVector vec
    vec = ch.ChVector[double](arr[0], arr[1], arr[2])
    return vec

cdef ch.ChQuaternion npArray_to_ChQuaternion(np.ndarray arr):
    cdef ch.ChQuaternion quat
    quat = ch.ChQuaternion[double](arr[0], arr[1], arr[2], arr[3])
    return quat


cdef class ChVector:
    """Cython class for ChVector
    (!) Does not use pointer
    """

    def __cinit__(self, double x, double y, double z):
        self.cppobj = ch.ChVector(x, y, z)
        
    cpdef double x(self):
        return self.cppobj.x()
    
    cpdef double y(self):
        return self.cppobj.y()
    
    cpdef double z(self):
        return self.cppobj.z()


cdef class ChQuaternion:
    """Cython class for ChQuaternion
    (!) Does not use pointer
    """

    def __cinit__(self, double e0, double e1, double e2, double e3):
        self.cppobj = ch.ChQuaternion(e0, e1, e2, e3)

    cpdef double e0(self):
        return self.cppobj.e0()
    
    cpdef double e1(self):
        return self.cppobj.e1()

    cpdef double e2(self):
        return self.cppobj.e2()

    cpdef double e3(self):
        return self.cppobj.e3()


cdef class ChCoordsys:
    """Cython class for ChCoordSys
    (!) Does not use pointer
    """

    def __cinit__(self, ChVector mv, ChQuaternion mq):
        self.cppobj = ch.ChCoordsys(mv.cppobj, mq.cppobj)

    cpdef np.ndarray getPos(self):
        return ChVector_to_npArray(self.cppobj.pos)

    cpdef np.ndarray getRot(self):
        return ChQuaternion_to_npArray(self.cppobj.rot)


cdef class ChFrame:
    """Cython class for ChFrame
    (!) Uses shared_ptr
    """

    def __cinit(self):
        if type(self) is ChFrame:
            self.sharedptr_chframe = make_shared[ch.ChFrame]()

    cpdef np.ndarray GetPos(self):
        return ChVector_to_npArray(deref(self.sharedptr_chframe).GetPos())

    cpdef void SetPos(self, ChVector mpos):
        deref(self.sharedptr_chframe).SetPos(mpos.cppobj)

    cpdef np.ndarray GetRot(self):
        return ChQuaternion_to_npArray(deref(self.sharedptr_chframe).GetRot())
    
    cpdef void SetRot(self, ChQuaternion rot):
        deref(self.sharedptr_chframe).SetRot(rot.cppobj)

    cpdef np.ndarray GetA(self):
        return ChMatrix33_to_npArray(deref(self.sharedptr_chframe).GetA())

    cpdef np.ndarray GetRotAxis(self):
        cdef ch.ChVector vec
        vec = deref(self.sharedptr_chframe).GetRotAxis()
        return ChVector_to_npArray(vec)

    cpdef double GetRotAngle(self):
        return deref(self.sharedptr_chframe).GetRotAngle()


cdef class ChFrameMoving(ChFrame):
    """Cython class for ChMovingFrame
    (!) Uses shared_ptr
    """

    def __cinit(self):
        if type(self) is ChFrameMoving:
            self.sharedptr_chframemoving = make_shared[ch.ChFrameMoving]()
            self.sharedptr_chframe = <shared_ptr[ch.ChFrame]> self.sharedptr_chframemoving

    cpdef np.ndarray GetPos_dt(self):
        return ChVector_to_npArray(deref(self.sharedptr_chframemoving).GetPos_dt())

    cpdef void SetPos_dt(self, ChVector mpos):
        deref(self.sharedptr_chframemoving).SetPos_dt(mpos.cppobj)

    cpdef np.ndarray GetPos_dtdt(self):
        return ChVector_to_npArray(deref(self.sharedptr_chframemoving).GetPos_dtdt())

    cpdef void SetPos_dtdt(self, ChVector mpos):
        deref(self.sharedptr_chframemoving).SetPos_dtdt(mpos.cppobj)

    cpdef np.ndarray GetRot_dt(self):
        return ChQuaternion_to_npArray(deref(self.sharedptr_chframemoving).GetRot_dt())

    cpdef void SetRot_dt(self, ChQuaternion mrot):
        deref(self.sharedptr_chframemoving).SetRot_dt(mrot.cppobj)

    cpdef np.ndarray GetRot_dtdt(self):
        return ChQuaternion_to_npArray(deref(self.sharedptr_chframemoving).GetRot_dtdt())

    cpdef void SetRot_dtdt(self, ChQuaternion mrot):
        deref(self.sharedptr_chframemoving).SetRot_dtdt(mrot.cppobj)

    cpdef np.ndarray GetWvel_loc(self):
        cdef ch.ChVector vec
        vec = deref(self.sharedptr_chframemoving).GetWvel_loc()
        return ChVector_to_npArray(vec)

    cpdef np.ndarray GetWacc_loc(self):
        cdef ch.ChVector vec
        vec = deref(self.sharedptr_chframemoving).GetWacc_loc()
        return ChVector_to_npArray(vec)


cdef class ChBodyFrame(ChFrameMoving):
    """Cython class for ChBodyFrame
    (!) Uses shared_ptr
    """

    def __cinit(self):
        # following commented does not work because abstract class
        pass
        # if type(self) is ChBodyFrame:
        #     self.sharedptr_chbodyframe = make_shared[ch.ChBodyFrame]()
        #     self.sharedptr_chframemoving = <shared_ptr[ch.ChFrameMoving]> self.sharedptr_chbodyframe
        #     self.sharedptr_chframe = <shared_ptr[ch.ChFrame]> self.sharedptr_chbodyframe


cdef class ChBody(ChBodyFrame):
    """Cython class for ChBody
    (!) Uses shared_ptr
    """

    def __cinit__(self):
        if type(self) is ChBody:
            self.sharedptr_chbody = make_shared[ch.ChBody]()
            self.sharedptr_chbodyframe = <shared_ptr[ch.ChBodyFrame]> self.sharedptr_chbody
            self.sharedptr_chframemoving = <shared_ptr[ch.ChFrameMoving]> self.sharedptr_chbody
            self.sharedptr_chframe = <shared_ptr[ch.ChFrame]> self.sharedptr_chbody

    cpdef void SetBodyFixed(self, bool state):
        deref(self.sharedptr_chbody).SetBodyFixed(state)

    cpdef void SetMaterialSurface(self, ChMaterialSurfaceSMC mat):
        deref(self.sharedptr_chbody).SetMaterialSurface(<shared_ptr[ch.ChMaterialSurface]> mat.sharedptr)


cdef class ChBodyEasyBox(ChBody):
    """Cython class for ChBodyEasyBox
    (!) Uses shared_ptr
    """

    def __cinit__(self, double Xsize, double Ysize, double Zsize, double mdensity, bool collide=True, bool visual_asset=False):
        if type(self) is ChBodyEasyBox:
            self.sharedptr_chbodyeasybox = make_shared[ch.ChBodyEasyBox](Xsize, Ysize, Zsize, mdensity, collide, visual_asset)
            self.sharedptr_chbody = <shared_ptr[ch.ChBody]> self.sharedptr_chbodyeasybox
            self.sharedptr_chbodyframe = <shared_ptr[ch.ChBodyFrame]> self.sharedptr_chbodyeasybox
            self.sharedptr_chframemoving = <shared_ptr[ch.ChFrameMoving]> self.sharedptr_chbodyeasybox
            self.sharedptr_chframe = <shared_ptr[ch.ChFrame]> self.sharedptr_chbodyeasybox
            # to remove below


cdef class ChMaterialSurfaceSMC:
    """cython class for chmaterialsurfacedem
    (!) uses shared_ptr
    """

    def __cinit__(self):
        self.sharedptr = make_shared[ch.ChMaterialSurfaceSMC]()

    cpdef void SetYoungModulus(self, float val):
        deref(self.sharedptr).SetYoungModulus(val)

    cpdef void SetPoissonRatio(self, float val):
        deref(self.sharedptr).SetPoissonRatio(val)

    cpdef void SetSfriction(self, float val):
        deref(self.sharedptr).SetSfriction(val)

    cpdef void SetKfriction(self, float val):
        deref(self.sharedptr).SetKfriction(val)

    cpdef void SetFriction(self, float val):
        deref(self.sharedptr).SetFriction(val)

    cpdef void SetRestitution(self, float val):
        deref(self.sharedptr).SetRestitution(val)

    cpdef void SetAdhesion(self, float val):
        deref(self.sharedptr).SetAdhesion(val)

    cpdef void SetAdhesionMultDMT(self, float val):
        deref(self.sharedptr).SetAdhesionMultDMT(val)

    cpdef void SetKn(self, float val):
        deref(self.sharedptr).SetKn(val)

    cpdef void SetKt(self, float val):
        deref(self.sharedptr).SetKt(val)

    cpdef void SetGn(self, float val):
        deref(self.sharedptr).SetGn(val)

    cpdef void SetGt(self, float val):
        deref(self.sharedptr).SetGt(val)

    cpdef double GetYoungModulus(self):
        return deref(self.sharedptr).GetYoungModulus()

    cpdef double GetPoissonRatio(self):
        return deref(self.sharedptr).GetPoissonRatio()

    cpdef double GetSfriction(self):
        return deref(self.sharedptr).GetSfriction()

    cpdef double GetKfriction(self):
        return deref(self.sharedptr).GetKfriction()

    cpdef double GetRestitution(self):
        return deref(self.sharedptr).GetRestitution()

    cpdef double GetAdhesion(self):
        return deref(self.sharedptr).GetAdhesion()

    cpdef double GetAdhesionMultDMT(self):
        return deref(self.sharedptr).GetAdhesionMultDMT()
    
    cpdef double GetKn(self):
        return deref(self.sharedptr).GetKn()

    cpdef double GetKt(self):
        return deref(self.sharedptr).GetKt()

    cpdef double GetGn(self):
        return deref(self.sharedptr).GetGn()

    cpdef double GetGt(self):
        return deref(self.sharedptr).GetGt()


cdef class ChContactSurfaceNodeCloud:
    """Cython class for ChContactSurfaceNodeCloud
    (!) Uses shared_ptr
    """

    def __cinit__(self):
        self.sharedptr = make_shared[ch.ChContactSurfaceNodeCloud]()

    cpdef void AddAllNodes(self, double point_radius):
        deref(self.sharedptr).AddAllNodes(point_radius)


# cdef class ChLinkPointPoint:
#     cdef shared_ptr[ch.ChLinkPointPoint] sharedptr
#     def __cinit__(self):
#         self.sharedptr = make_shared[ch.ChLinkPointPoint]()
#     def Initialize(shared_ptr[ChNodeFEAxyz] anodeA, shared_ptr[ChNodeFEAxyz])
#         deref(self.sharedptr).Initialize(<shared_ptr[ch.ChNodeFEAxyz]> anodeA, <shared_ptr[ch.ChNodeFEAxyz]> anodeB)





cdef extern from "swigpyobject.h":
    ctypedef struct SwigPyObject:
        void *ptr



cdef class ChBodyAddedMass:
    """Cython class for ChBodyAddedMass
    (!) Uses shared_ptr
    """

    def __cinit__(self):
        self.sharedptr = make_shared[ch.ChBodyAddedMass]()
        self.sharedptr_chbody = <shared_ptr[ch.ChBody]> self.sharedptr
        self.thisptr = self.sharedptr.get()
        self.ChBodySWIG = chrono.ChBody()
        cdef SwigPyObject *swig_obj = <SwigPyObject*>self.ChBodySWIG.this
        self.bodyptr = <ch.ChBody*?> &self.thisptr
        swig_obj.ptr = self.bodyptr
        # print(self.ChBodyy.GetPos())

        # cdef ch.ChVe *mycpp_ptr = <ch.ChVector*?>swig_obj.ptr
        # cdef ch.ChVector my_instance = deref(mycpp_ptr)

    cdef void SetMfullmass(self, ch.ChMatrixDynamic Mfullmass_in):
        self.thisptr.SetMfullmass(Mfullmass_in)

    cdef void SetInvMfullmass(self, ch.ChMatrixDynamic inv_Mfullmass_in):
        self.thisptr.SetInvMfullmass(inv_Mfullmass_in)

    cpdef np.ndarray GetInertia(self):
        return ConstChMatrix33_to_npArray(self.thisptr.GetInertia())

    cpdef void SetInertiaXX(self, ChVector iner):
        self.thisptr.SetInertiaXX(iner.cppobj)

    cpdef void SetInertiaXY(self, ChVector iner):
        self.thisptr.SetInertiaXY(iner.cppobj)

    cpdef void SetMass(self, double newmass):
        self.thisptr.SetMass(newmass)

    cpdef double GetMass(self):
        return self.thisptr.GetMass()

    cpdef double GetPos(self):
        return self.thisptr.GetPos().x()
