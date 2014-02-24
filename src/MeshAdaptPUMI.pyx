# distutils: language = c++

from cpython.ref cimport PyObject
import numpy as np
cimport numpy as np

cdef extern from "mesh.h":
    struct Mesh:
       pass

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
        int UpdateMaterialArrays(Mesh&, int, int)
        int TransferSolutionToPUMI(double*, int)

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
    def UpdateMaterialArrays(self, cmesh, bdryId, geomTag):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.UpdateMaterialArrays(cmesh_ptr.mesh, bdryId, geomTag)
    def TransferSolutionToPUMI(self, np.ndarray[double,ndim=1,mode="c"] input not None, int index)
        return self.thisptr.TransferSolutionToPUMI(&input[0], index)
