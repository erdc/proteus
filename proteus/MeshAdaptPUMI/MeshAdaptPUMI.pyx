# distutils: language = c++

from cpython.ref cimport PyObject
cimport numpy as np
import numpy as np
from ..Profiling import logEvent
from libcpp.string cimport string

cdef extern from "mesh.h":
    struct Mesh:
       pass

cdef extern from "cmeshToolsModule.h":
    ctypedef struct CMesh:
        Mesh mesh


cdef extern from "MeshAdaptPUMI/MeshAdaptPUMI.h":
    cdef cppclass MeshAdaptPUMIDrvr:
        MeshAdaptPUMIDrvr(double, double, int, char*, char*,char*,double,double)
        int numIter, numAdaptSteps
        string size_field_config, adapt_type_config
        bint isReconstructed
        int loadModelAndMesh(char *, char*)
        int getSimmetrixBC()
        int reconstructFromProteus(Mesh&,Mesh&,int)
        int constructFromSerialPUMIMesh(Mesh&)
        int constructFromParallelPUMIMesh(Mesh&, Mesh&)
        int updateMaterialArrays(Mesh&,int, int, int)
        int updateMaterialArrays(Mesh&)
        int transferFieldToPUMI(char*, double*, int, int)
        int transferFieldToProteus(char*, double*, int, int)
        int transferPropertiesToPUMI(double*, double*,double*)
        int transferModelInfo(int*,int*,int*,int*,int*,int*,int)
        int transferBCtagsToProteus(int*, int, int*, int*,double*)
        int transferBCsToProteus()
        int adaptPUMIMesh()
        int dumpMesh(Mesh&)
        int willAdapt()
        int getERMSizeField(double)
        double getMinimumQuality()
        double getTotalMass()
        double getMPvalue(double,double, double)
        void get_local_error(double)

cdef class MeshAdaptPUMI:
    cdef MeshAdaptPUMIDrvr *thisptr
    cdef double hmax, hmin
    cdef int numIter, numAdaptSteps
    cdef bint isReconstructed
    def __cinit__(self, hmax=100.0, hmin=1e-8, numIter=10, sfConfig="ERM",maType="isotropic",logType="off",targetError=0,targetElementCount=0):
        logEvent("MeshAdaptPUMI: hmax = {0} hmin = {1} numIter = {2}".format(hmax,hmin,numIter))
        self.thisptr = new MeshAdaptPUMIDrvr(hmax, hmin, numIter, sfConfig,maType,logType,targetError,targetElementCount)
    def __dealloc__(self):
        del self.thisptr
    def size_field_config(self):
        return self.thisptr.size_field_config
    def adapt_type_config(self):
        return self.thisptr.adapt_type_config
    def numIter(self):
        return self.thisptr.numIter
    def isReconstructed(self):
        return self.thisptr.isReconstructed
    def loadModelAndMesh(self, geomName, meshName):
        return self.thisptr.loadModelAndMesh(geomName, meshName)
    def reconstructFromProteus(self,cmesh,global_cmesh,hasModel=0):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        cdef CMesh* global_cmesh_ptr = <CMesh*>global_cmesh
        return self.thisptr.reconstructFromProteus(cmesh_ptr.mesh,global_cmesh_ptr.mesh,hasModel)
    def constructFromSerialPUMIMesh(self, cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.constructFromSerialPUMIMesh(cmesh_ptr.mesh)
    def constructFromParallelPUMIMesh(self, cmesh, subdomain_cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        cdef CMesh* subdomain_cmesh_ptr = <CMesh*>subdomain_cmesh
        return self.thisptr.constructFromParallelPUMIMesh(cmesh_ptr.mesh, subdomain_cmesh_ptr.mesh)
    def updateMaterialArrays(self, cmesh, dim=None,bdryId=None, geomTag=None):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        if(dim==None):
            return self.thisptr.updateMaterialArrays(cmesh_ptr.mesh)
        else:
            return self.thisptr.updateMaterialArrays(cmesh_ptr.mesh,dim, bdryId, geomTag)
    def transferFieldToPUMI(self, name, np.ndarray[np.double_t,ndim=2,mode="c"] inArray):
        inArray = np.ascontiguousarray(inArray)
        return self.thisptr.transferFieldToPUMI(name, &inArray[0,0], inArray.shape[1], inArray.shape[0])
    def transferFieldToProteus(self, name, np.ndarray[np.double_t,ndim=2,mode="c"] outArray):
        outArray = np.ascontiguousarray(outArray)
        return self.thisptr.transferFieldToProteus(name, &outArray[0,0], outArray.shape[1], outArray.shape[0])
    def transferPropertiesToPUMI(self, np.ndarray[np.double_t,ndim=1,mode="c"] rho, np.ndarray[np.double_t,ndim=1,mode="c"] nu, np.ndarray[np.double_t,ndim=1,mode="c"] g):
        rho = np.ascontiguousarray(rho)
        nu = np.ascontiguousarray(nu)
        g = np.ascontiguousarray(g)
        return self.thisptr.transferPropertiesToPUMI(&rho[0],&nu[0],&g[0])
    def transferModelInfo(self, np.ndarray[int,ndim=1,mode="c"] numModelEntities,
                                np.ndarray[int,ndim=2,mode="c"] edges,
                                np.ndarray[int,ndim=2,mode="c"] faces,
                                np.ndarray[int,ndim=2,mode="c"] meshVertex2Model,
                                np.ndarray[int,ndim=2,mode="c"] meshEdge2Model,
                                np.ndarray[int,ndim=2,mode="c"] meshBoundary2Model
    ):
        numModelEntities = np.ascontiguousarray(numModelEntities)
        edges = np.ascontiguousarray(edges)
        nMaxSegments = faces.shape[1]
        faces = np.ascontiguousarray(faces)
        meshVertex2Model = np.ascontiguousarray(meshVertex2Model)
        meshEdge2Model = np.ascontiguousarray(meshEdge2Model)
        meshBoundary2Model = np.ascontiguousarray(meshBoundary2Model)
        return self.thisptr.transferModelInfo(<int*> numModelEntities.data,<int*> edges.data,<int *> faces.data,<int*> meshVertex2Model.data,<int*> meshEdge2Model.data,<int*> meshBoundary2Model.data,nMaxSegments)
    #def transferBCtagsToProteus(self, np.ndarray[int,ndim=2,mode="c"] tagArray, int idx, np.ndarray[int,ndim=1,mode="c"] ebN, np.ndarray[int, ndim=2, mode="c"] eN_global, np.ndarray[np.double_t,ndim=2,mode="c"] fluxBC):
    #    tagArray = np.ascontiguousarray(tagArray)
    #    ebN = np.ascontiguousarray(ebN)
    #    eN_global = np.ascontiguousarray(eN_global)
    #    fluxBC = np.ascontiguousarray(fluxBC)
    #    return self.thisptr.transferBCtagsToProteus(&tagArray[0,0],idx,&ebN[0],&eN_global[0,0],&fluxBC[0,0])
    #def transferBCsToProteus(self):
    #    return self.thisptr.transferBCsToProteus()
    def adaptPUMIMesh(self):
        return self.thisptr.adaptPUMIMesh()
    def dumpMesh(self, cmesh):
        cdef CMesh* cmesh_ptr = <CMesh*>cmesh
        return self.thisptr.dumpMesh(cmesh_ptr.mesh)
    def getERMSizeField(self, err_total):
        return self.thisptr.getERMSizeField(err_total);
    def getMPvalue(self,field_val,val_0,val_1):
        return self.thisptr.getMPvalue(field_val,val_0,val_1)
    def willAdapt(self):
        return self.thisptr.willAdapt()
    def get_local_error(self):
        errTotal=0.0;
        self.thisptr.get_local_error(errTotal)
        return errTotal
