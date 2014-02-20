# distutils: language = c++

from cpython.ref cimport PyObject

cdef extern from "mesh.h":
    struct Mesh:
       pass
#       int *elementOffsets_subdomain_owned
#       int *elementBoundaryOffsets_subdomain_owned
#       int *edgeOffsets_subdomain_owned
#       int *nodeOffsets_subdomain_owned
#       int *elementNumbering_subdomain2global
#       int *elementBoundaryNumbering_subdomain2global
#       int *edgeNumbering_subdomain2global
#       int *nodeNumbering_subdomain2global          

cdef extern from "cmeshToolsModule.h":
    ctypedef struct CMesh:
        Mesh mesh


cdef extern from "MeshAdaptPUMI/MeshAdaptPUMI.h":
    cdef cppclass MeshAdaptPUMIDrvr:
        MeshAdaptPUMIDrvr()
        int helloworld(char *)
        int readGeomModel(char *)
        int readPUMIMesh(char *)
        int initProteusMesh(Mesh&)
        int ConstructFromSerialPUMIMesh(Mesh&)
        int ConstructFromParallelPUMIMesh(Mesh&, Mesh&)


cdef class MeshAdaptPUMI:
    cdef MeshAdaptPUMIDrvr *thisptr
    def __cinit__(self):
        self.thisptr = new MeshAdaptPUMIDrvr()
    def __dealloc__(self):
        del self.thisptr
    def helloworld(self, string):
        return self.thisptr.helloworld(string)
    def readGeomModel(self, geomName):
        return self.thisptr.readGeomModel(geomName)
    def readPUMIMesh(self, meshName):
        return self.thisptr.readPUMIMesh(meshName)
    def initProteusMesh(self, cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.initProteusMesh(cmesh_ptr.mesh)
    def ConstructFromSerialPUMIMesh(self, cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.ConstructFromSerialPUMIMesh(cmesh_ptr.mesh)
    def ConstructFromParallelPUMIMesh(self, cmesh, subdomain_cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        cdef CMesh* subdomain_cmesh_ptr = <CMesh*>subdomain_cmesh
        return self.thisptr.ConstructFromParallelPUMIMesh(cmesh_ptr.mesh, subdomain_cmesh_ptr.mesh)
      
#        return cmesh_ptr.mesh.elementOffsets_subdomain_owned, cmesh_ptr.mesh.elementNumbering_subdomain2global, cmesh_ptr.mesh.nodeOffsets_subdomain_owned, cmesh_ptr.mesh.nodeNumbering_subdomain2global, cmesh_ptr.mesh.elementBoundaryOffsets_subdomain_owned, cmesh_ptr.mesh.elementBoundaryNumbering_subdomain2global, cmesh_ptr.mesh.edgeOffsets_subdomain_owned, cmesh_ptr.mesh.edgeNumbering_subdomain2globali




