# distutils: language = c++

from cpython.ref cimport PyObject
cimport numpy as np
import numpy as np
from ..Profiling import logEvent

cdef extern from "mesh.h":
    struct Mesh:
       pass

cdef extern from "cmeshToolsModule.h":
    ctypedef struct CMesh:
        Mesh mesh


cdef extern from "MeshAdaptPUMI/MeshAdaptPUMI.h":
    cdef cppclass MeshAdaptPUMIDrvr:
        MeshAdaptPUMIDrvr(double, double, int, char*)
        int numIter, numAdaptSteps
        int loadModelAndMesh(char *, char*)
        int constructFromSerialPUMIMesh(Mesh&)
        int constructFromParallelPUMIMesh(Mesh&, Mesh&)
        int updateMaterialArrays(Mesh&, int, int)
        int transferSolutionToPUMI(double*, int, int)
        int transferSolutionToProteus(double*, int, int)
        int transferFieldToPUMI(char*, double*, int, int)
        int transferPropertiesToPUMI(double*, double*)
        int transferBCtagsToProteus(int*, int, int*, int*,double*)
        int transferBCsToProteus()
        int adaptPUMIMesh()
        int dumpMesh(Mesh&)
        int getERMSizeField(double)
        double getMinimumQuality()

cdef class MeshAdaptPUMI:
    cdef MeshAdaptPUMIDrvr *thisptr
    cdef double hmax, hmin
    cdef int numIter, numAdaptSteps
    def __cinit__(self, hmax=100.0, hmin=1e-8, numIter=10, sfConfig="farhad"):
        logEvent("MeshAdaptPUMI: hmax = {0} hmin = {1} numIter = {2}".format(hmax,hmin,numIter))
        self.thisptr = new MeshAdaptPUMIDrvr(hmax, hmin, numIter, sfConfig)
    def __dealloc__(self):
        del self.thisptr
    def loadModelAndMesh(self, geomName, meshName):
        return self.thisptr.loadModelAndMesh(geomName, meshName)
    def constructFromSerialPUMIMesh(self, cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.constructFromSerialPUMIMesh(cmesh_ptr.mesh)
    def constructFromParallelPUMIMesh(self, cmesh, subdomain_cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        cdef CMesh* subdomain_cmesh_ptr = <CMesh*>subdomain_cmesh
        return self.thisptr.constructFromParallelPUMIMesh(cmesh_ptr.mesh, subdomain_cmesh_ptr.mesh)
    def updateMaterialArrays(self, cmesh, bdryId, geomTag):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.updateMaterialArrays(cmesh_ptr.mesh, bdryId, geomTag)
    def transferSolutionToPUMI(self, np.ndarray[np.double_t,ndim=2,mode="c"] inArray):
        inArray = np.ascontiguousarray(inArray)
        return self.thisptr.transferSolutionToPUMI(&inArray[0,0], inArray.shape[0], inArray.shape[1])
    def transferSolutionToProteus(self, np.ndarray[np.double_t,ndim=2,mode="c"] outArray):
        outArray = np.ascontiguousarray(outArray)
        return self.thisptr.transferSolutionToProteus(&outArray[0,0], outArray.shape[0], outArray.shape[1])
    def transferFieldToPUMI(self, name, np.ndarray[np.double_t,ndim=2,mode="c"] inArray):
        inArray = np.ascontiguousarray(inArray)
        return self.thisptr.transferFieldToPUMI(name, &inArray[0,0], inArray.shape[1], inArray.shape[0])
    def transferPropertiesToPUMI(self, np.ndarray[np.double_t,ndim=1,mode="c"] rho, np.ndarray[np.double_t,ndim=1,mode="c"] nu):
        rho = np.ascontiguousarray(rho)
        nu = np.ascontiguousarray(nu)
        return self.thisptr.transferPropertiesToPUMI(&rho[0],&nu[0])
    def transferBCtagsToProteus(self, np.ndarray[int,ndim=2,mode="c"] tagArray, int idx, np.ndarray[int,ndim=1,mode="c"] ebN, np.ndarray[int, ndim=2, mode="c"] eN_global, np.ndarray[np.double_t,ndim=2,mode="c"] fluxBC):
        tagArray = np.ascontiguousarray(tagArray)
        ebN = np.ascontiguousarray(ebN)
        eN_global = np.ascontiguousarray(eN_global)
        fluxBC = np.ascontiguousarray(fluxBC)
        return self.thisptr.transferBCtagsToProteus(&tagArray[0,0],idx,&ebN[0],&eN_global[0,0],&fluxBC[0,0])
    def transferBCsToProteus(self):
        return self.thisptr.transferBCsToProteus()
    def adaptPUMIMesh(self):
        return self.thisptr.adaptPUMIMesh()
    def dumpMesh(self, cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.dumpMesh(cmesh_ptr.mesh)
    def getERMSizeField(self, err_total):
        return self.thisptr.getERMSizeField(err_total);
