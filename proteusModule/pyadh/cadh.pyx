cimport cadh
import os
import numpy
cimport numpy
ctypedef numpy.double_t DTYPE_t
#numpy.intc_t not in cython's numpy.pxd 
ctypedef int ITYPE_t
import cmeshTools
import MeshTools
#for creating numpy arrays from adh data pointers
cdef extern from "Python.h":
    cdef void Py_INCREF(object)
    ctypedef struct PyTypeObject:
        pass
cdef extern from "numpy/arrayobject.h":
    cdef void import_array()
    cdef numpy.ndarray PyArray_SimpleNewFromData(int nd, numpy.npy_intp* dims, int typenum, void* data)
    #
    PyTypeObject PyArray_Type
    ctypedef struct PyArray_Descr:
        pass
    PyArray_Descr * PyArray_DescrNewFromType(int type_num)
    cdef numpy.ndarray PyArray_NewFromDescr(PyTypeObject *subtype,
                                            numpy.dtype, 
                                            int nd,
                                            numpy.npy_intp* dims,
                                            numpy.npy_intp* strides,
                                            void* data,
                                            int flags,
                                            object parent)

    #set flags here to force write access?
    int NPY_WRITEABLE
    int NPY_CARRAY
cdef extern from "share_extern.h":
    pass

cdef class cADH_InputTranslator:
    cdef cadh.ADH_Input c_adh
    def __cinit__(self,filebase):
        cadh.readBC(&self.c_adh,filebase)

    
cdef class cADH_NumericalSolution:
    cdef cadh.ADH_NumericalSolution c_adh
    cdef object communicator
    def __cinit__(self,filename,runname,comm):
        cdef int argc = 3
        cdef char* argv[3]
        cdef int failure = 0

        argv[0] = "adh"
        argv[1] = <char*>filename
        argv[2] = <char*>runname

        pyadh_packages = os.getenv('PYADH_PACKAGES',os.getenv('HOME')+'/src/pyadh-packages')
        preadh = pyadh_packages+'/adh/bin/pre_adh'
        command = preadh+' '+filename
        print "cadh calling pre_adh as %s" % command
        os.system(command)
        
        failure = cadh.ADH_NumericalSolution_init(&self.c_adh,argc,argv)
        self.communicator = comm
#    def __init__(self,filename,runname,comm):
#        self.communicator = comm
#     property comm:
#         def __get__(self):
#             return self.communicator
#         def __set__(self,value):
#             self.communicator = value
#         def __del__(self):
#             pass
    def __dealloc__(self):
        cdef int failure = 0
        print "calling cADH_NumericalSolution_dealloc from pyx"
        failure = cadh.ADH_NumericalSolution_dealloc(&self.c_adh)
        print "cADH_NumericalSolution_dealloc from pyx returned %s" % failure

    def step(self):
        cdef int failure = 0
        failure = cadh.ADH_NumericalSolution_step(&self.c_adh)
        return failure
    def stepTaken(self):
        cdef int failure = 0
        failure = cadh.ADH_NumericalSolution_stepTaken(&self.c_adh)
        return failure

    def calculateSolution(self):
        cdef int failure = 0
        failure = cadh.ADH_NumericalSolution_calculateSolution(&self.c_adh)
        return self

cdef class cADH_OneLevelTransport:
    cdef cadh.ADH_OneLevelTransport transport
    cdef object numericalSolution
    def __cinit__(self,numericalSolution):#,comm):
        cadh.ADH_OneLevelTransport_init(&self.transport)
        self.numericalSolution = numericalSolution
    def __dealloc__(self):
        cdef int failure = 0
        print "calling cadh.ADH_OneLevelTransport_dealloc"
        failure = cadh.ADH_OneLevelTransport_dealloc(&self.transport)
        print "back from cadh.ADH_OneLevelTransport_dealloc failure= %s" % failure
    def realloc(self):
        cadh.realloc_fe_main(&self.transport)
    def updateStorage(self):
        cadh.ADH_OneLevelTransport_updateStorage(&self.transport)
    def solve(self):
        cdef int failure = 0
        failure = cadh.ADH_OneLevelTransport_solve(&self.transport)
        return failure
    def getSolutionComponentDimension(self, int ci):
        """
        size of component ci
        That is it is 1 for a scalar variable like pressure
        and 2 for the 2D shallow water velocitiy 
        """
        return cadh.get_ADH_solutionComponentDimension(ci)
    def getFirstSolutionComponentForVectorComponent(self, int ci):
        return cadh.get_ADH_firstSolutionComponentForVectorComponent(ci)
    def getVectorSolutionComponentIds(self):
        """
        return which  components are logically part of the model's vector unknown0
        """
        cdef int nc  = cadh.get_ADH_NumberOfComponents()
        cdef list vectorComponents = []
        cdef int ndim
        for ci in range(nc):
            ndim = cadh.get_ADH_solutionComponentDimension(ci)
            if ndim > 1:
                vectorComponents.append(ci)
        return vectorComponents
    def getSolutionComponentView(self,int ci):
        """
        get a strided view of a vector solution
        """
        cdef int nd = 1
        cdef int flags = NPY_WRITEABLE
        cdef numpy.npy_intp dims[1]
        cdef numpy.npy_intp byte_strides[1]
        cdef numpy.dtype dt = numpy.dtype('d')
        #need a solution dimension function
        cdef int ndof_base = cadh.get_ADH_nNodes_global()
        cdef int nc  = cadh.get_ADH_NumberOfComponents()
        #can be > 1 if represents a vector valued quantity
        #eg shallow water components (0,1) are (u,v)
        # si this would return 2 if ci = 0 or 1
        cdef int ndim  = cadh.get_ADH_solutionComponentDimension(ci)
        assert ndim > 0
        cdef int ci_base = cadh.get_ADH_firstSolutionComponentForVectorComponent(ci)
        cdef int ci_adh = ci-ci_base
        #get glocal solution from adh
        cdef void * u = cadh.get_ADH_solutionComponent(ci)
        assert u != NULL
        #now return a view of u that is only a single component 
        dims[0] = ndof_base
        byte_strides[0] = ndim*sizeof(double)
        cdef numpy.ndarray pysol
        Py_INCREF(dt)
        if (ndim > 1):
            assert ci_adh >= 0
            u += ci_adh*sizeof(double)
            pysol = PyArray_NewFromDescr(&PyArray_Type,
                                         dt,
                                         nd,
                                         dims,
                                         byte_strides,
                                         u,
                                         flags,
                                         <object>NULL)
            return pysol
        else:
            pysol = PyArray_NewFromDescr(&PyArray_Type,
                                         dt,
                                         nd,
                                         dims,
                                         NULL,
                                         u,
                                         flags,
                                         <object>NULL)
            return pysol
            
    def getSolutionComponent(self,int ci):
        """
        return component ci from adh solution
        if component ci is part of a vector valued quantity,
          then the entire vector valued component is returned

        For example, if the shallow water
          componets are (u,v,h) and the calling routine passes in
          ci=0 or ci=1, it will get \vec v= (u,v) while ci=2 will return h
        """
        cdef int nd = 1
        cdef numpy.npy_intp dims[1]
        #need a solution dimension function
        cdef int ndof_base = cadh.get_ADH_nNodes_global()
        cdef int ndim  = cadh.get_ADH_solutionComponentDimension(ci)
        cdef int flags = NPY_WRITEABLE
        assert ndim > 0;
        cdef int ndof  = ndim*ndof_base
        #get glocal solution from adh
        cdef void * u = cadh.get_ADH_solutionComponent(ci)
        assert u != NULL
        cdef numpy.dtype dt
        if ndim == 1:
            dt = numpy.dtype('d')
        else:
            #numpy rec array
            format = [('%s' % i,numpy.double) for i in range(ndim)]
            dt = numpy.dtype(format)
        Py_INCREF(dt)
        dims[0] = ndof_base
        print "in cadh getSolutionComponent(%s) dims[0]= %s dt= %s" % (ci,dims[0],dt)
        cdef numpy.ndarray pysol = PyArray_NewFromDescr(&PyArray_Type,
                                                        dt,
                                                        nd,
                                                        dims,
                                                        NULL,
                                                        <void*>u,
                                                        flags,
                                                        <object>NULL)
        
        return pysol

    def meshAdaptedForNextStep(self):
        """

        """
        cdef int adapted = cadh.get_ADH_ADAPTION_flag()
        return adapted
    property nc:
        """
        number of components in system
        """
        def __get__(self):
            return cadh.get_ADH_NumberOfComponents()
        def __set__(self,value):
            assert False, "can't set adh number of componets yet"
        def __del__(self):
            pass
    property nDOF_global:
        """
        total number of degrees of freedom
        """
        def __get__(self):
            return cadh.get_ADH_NumberOfDegreesOfFreedom()
        def __set__(self,value):
            assert False, "can't set adh number of dof yet"
        def __del__(self):
            pass
    property t:
        """
        new time level 
        """
        def __get__(self):
            return cadh.get_ADH_newTimeLevel()
        def __set__(self,value):
            assert False, "can't set adh time level yet"
        def __del__(self):
            pass
    property t0:
        """
        starting time
        """
        def __get__(self):
            return cadh.get_ADH_tinitial()
        def __set__(self,value):
            assert False, "can't set adh initial time level yet"
        def __del__(self):
            pass
        
    property du:
        """
        Newton increment
        """
        def __get__(self):
            cdef int nd = 1
            cdef numpy.npy_intp dims[1]
            cdef int ndof_base = cadh.get_ADH_nNodes_global()
            cdef int nc  = cadh.get_ADH_NumberOfComponents()
            cdef int ndof  = nc*ndof_base
            dims[0] = ndof
            pydu = PyArray_SimpleNewFromData(nd,dims,numpy.NPY_DOUBLE,<void*>self.transport.sol)
            return pydu
        def __set__(self,numpy.ndarray[DTYPE_t,ndim=1] value):
            cdef int i
            cdef int ndof_base = cadh.get_ADH_nNodes_global()
            cdef int nc  = cadh.get_ADH_NumberOfComponents()
            cdef int ndof  = nc*ndof_base
            assert ndof == value.shape[0]
            for i in range(ndof):
                self.transport.sol[i] = value[i]
        def __del__(self):
            pass
    property residual:
        """
        nonlinear residual
        """
        def __get__(self):
            cdef int nd = 1
            cdef numpy.npy_intp dims[1]
            cdef int ndof_base = cadh.get_ADH_nNodes_global()
            cdef int nc  = cadh.get_ADH_NumberOfComponents()
            cdef int ndof  = nc*ndof_base
            dims[0] = ndof
            pyr = PyArray_SimpleNewFromData(nd,dims,numpy.NPY_DOUBLE,<void*>self.transport.residual)
            return pyr
        def __set__(self,numpy.ndarray[DTYPE_t,ndim=1] value):
            cdef int i
            cdef int ndof_base = cadh.get_ADH_nNodes_global()
            cdef int nc  = cadh.get_ADH_NumberOfComponents()
            cdef int ndof  = nc*ndof_base
            assert ndof == value.shape[0]
            for i in range(ndof):
                self.transport.residual[i] = value[i]
        def __del__(self):
            pass

cdef class cADH_Mesh:
    cdef int c_nElements_global
    cdef int c_nNodes_global
    cdef int c_nSpace_global
    cdef int c_nNodes_element
    #
    
    def __cinit__(self):
        self.c_nElements_global = cadh.get_ADH_nElements_global()
        self.c_nNodes_global    = cadh.get_ADH_nNodes_global()
        self.c_nSpace_global    = cadh.get_ADH_nSpace_global()
        self.c_nNodes_element   = cadh.get_ADH_nNodes_element()

    property nElements_global:
        """
        total number of elements in domain
        """
        def __get__(self):
            self.c_nElements_global = cadh.get_ADH_nElements_global()
            return self.c_nElements_global
        def __set__(self,value):
            assert False, "can't set adh nElements_global yet"
        def __del__(self):
            pass
    property nNodes_global:
        def __get__(self):
            self.c_nNodes_global = cadh.get_ADH_nNodes_global()
            return self.c_nNodes_global
        def __set__(self,value):
            assert False, "can't set adh nNodes_global yet"
        def __del__(self):
            pass
    property nSpace_global:
        def __get__(self):
            self.c_nSpace_global = cadh.get_ADH_nSpace_global()
            return self.c_nSpace_global
        def __set__(self,value):
            assert False, "can't set adh nSpace_global yet"
        def __del__(self):
            pass
    property nNodes_element:
        def __get__(self):
            self.c_nNodes_element = cadh.get_ADH_nNodes_element()
            return self.c_nNodes_element
        def __set__(self,value):
            assert False, "can't set adh nSpace_global yet"
        def __del__(self):
            pass
    ##todo add properties that grab shallow copies of arrays
    def update(self):
        """
        update mesh info in case mesh has changed
        """
        self.c_nElements_global = cadh.get_ADH_nElements_global()
        self.c_nNodes_global    = cadh.get_ADH_nNodes_global()
        self.c_nSpace_global    = cadh.get_ADH_nSpace_global()
        self.c_nNodes_element   = cadh.get_ADH_nNodes_element()
        
    #
    def get_nodeArray(self, numpy.ndarray[DTYPE_t,ndim=2] nodeArray):
        cadh.get_nodeArray(<double*>nodeArray.data)
        #mwf debug
        #cdef int nN,I
        #for nN in range(self.c_nNodes_global):
        #    for I in range(self.c_nSpace_global):
        #        print "out nodeArray[%s,%s]= %s" % (nN,I,nodeArray[nN,I])
    def get_elementNodesArray(self, numpy.ndarray[ITYPE_t,ndim=2] elementNodesArray):
        cadh.get_elementNodesArray(<int*>elementNodesArray.data)
        #mwf debug
        #cdef int nE,nN
        #for nE in range(self.c_nElements_global):
        #    for nN in range(self.c_nNodes_element):
        #        print "out elementNodesArray[%s,%s]= %s" % (nE,nN,elementNodesArray[nE,nN])
    
    def get_elementMaterialTypes(self, numpy.ndarray[ITYPE_t,ndim=1] elementMaterialTypes):
        cadh.get_elementMaterialTypes(<int*>elementMaterialTypes.data)
    def get_nodeMaterialTypes(self, numpy.ndarray[ITYPE_t,ndim=1] nodeMaterialTypes):
        cadh.get_nodeMaterialTypes(<int*>nodeMaterialTypes.data)
    
    def generateMeshToolsMesh(self):
        print "in cadh generateMeshToolsMesh"
        if self.c_nSpace_global == 2:
            mesh = MeshTools.TriangularMesh()
        elif self.c_nSpace_global == 3:
            mesh = MeshTools.TetrahedralMesh()
        else:
            mesh = MeshTools.EdgeMesh()
        mesh.cmesh = cmeshTools.CMesh()
        
        cmeshTools.allocateNodeAndElementNodeDataStructures(mesh.cmesh,
                                                            self.c_nElements_global,
                                                            self.c_nNodes_global,
                                                            self.c_nNodes_element)
        
        (mesh.nElements_global,mesh.nNodes_global,mesh.nNodes_element,
         mesh.elementNodesArray,mesh.elementMaterialTypes,
         mesh.nodeMaterialTypes,mesh.nodeArray)  = cmeshTools.buildLevel0PythonMeshInterface(mesh.cmesh)
        
        self.get_nodeArray(mesh.nodeArray)
        self.get_elementNodesArray(mesh.elementNodesArray)
        self.get_elementMaterialTypes(mesh.elementMaterialTypes)
        self.get_nodeMaterialTypes(mesh.nodeMaterialTypes)

        cmeshTools.constructElementBoundaryElementsArray(mesh.cmesh)
        if self.c_nSpace_global == 2:
            cmeshTools.allocateGeometricInfo_triangle(mesh.cmesh)
            cmeshTools.computeGeometricInfo_triangle(mesh.cmesh)
        elif self.c_nSpace_global == 3:
            cmeshTools.allocateGeometricInfo_tetrahedron(mesh.cmesh)
            cmeshTools.computeGeometricInfo_tetrahedron(mesh.cmesh)
        else:
            cmeshTools.allocateGeometricInfo_edge(mesh.cmesh)
            cmeshTools.computeGeometricInfo_edge(mesh.cmesh)

        mesh.buildFromC(mesh.cmesh)
        
        return mesh

cdef class cADH_PETSc_Interface:
    cdef cadh.ADH_PETSc_Interface cadh_petsc_interface
    cdef cADH_OneLevelTransport oneLevelTransport
    def __cinit__(self,oneLevelTransport):
        self.oneLevelTransport = <cADH_OneLevelTransport?>oneLevelTransport
        failure = cadh.ADH_PETSc_Interface_init(&self.cadh_petsc_interface,&self.oneLevelTransport.transport)
    def __dealloc__(self):
        failure = cadh.ADH_PETSc_Interface_dealloc(&self.cadh_petsc_interface)
    def update(self):
        failure = cadh.ADH_PETSc_Interface_update(&self.cadh_petsc_interface,&self.oneLevelTransport.transport)

#necessary for Numpy C api calls
import_array()
