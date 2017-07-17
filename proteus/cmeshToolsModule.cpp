#include "cmeshToolsModule.h"
#include "mesh.h"
#include <algorithm>
/* //mwf added for python wrapper structures */
/* #include "flcbdfWrappersModule.h" */
#include <valarray>
#include <iostream>
#include <map>
#include <set>
/** \file cmeshToolsModule.cpp
    \defgroup cmeshTools cmeshTools
    \brief Python interface to mesh library 
    @{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))


//mwf
// typedef struct
// {
//   PyObject_HEAD
//   Mesh mesh;
// } CMesh;

// #define MESH(p) ((CMesh*)p)->mesh

extern "C"
{
static PyObject*
CMesh_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  CMesh *self;
  self = (CMesh *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
CMesh_init(CMesh *self, PyObject *args, PyObject *kwds)
{
  initializeMesh(self->mesh);
  return 0;
}

static  void
CMesh_dealloc(CMesh* self)
{
  deleteMesh(self->mesh);
}

static PyTypeObject CMeshType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "cmeshTools.CMesh",             /*tp_name*/
  sizeof(CMesh), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)CMesh_dealloc,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,        /*tp_flags*/
  "CMesh objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  0,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)CMesh_init,      /* tp_init */
  0,                         /* tp_alloc */
  CMesh_new,                 /* tp_new */
};

PyObject* CMesh_FromMesh(Mesh* mesh)
{
  CMesh* self;
  self = PyObject_NEW(CMesh, &CMeshType);
  self->mesh= *mesh;
  return (PyObject*)self;
}

static PyObject* cmeshToolsBuildPythonMeshInterface(PyObject* self,
                                                    PyObject* args)
{
  PyObject *cmesh,
    *elementNodesArray,
    *nodeElementsArray,
    *nodeElementOffsets,
    *elementNeighborsArray,
    *elementBoundariesArray,
    *elementBoundaryNodesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *interiorElementBoundariesArray,
    *exteriorElementBoundariesArray,
    *edgeNodesArray,
    *nodeStarArray,
    *nodeStarOffsets,
    *elementMaterialTypes,
    *elementBoundaryMaterialTypes,
    *nodeMaterialTypes,
    *nodeArray,
    *elementDiametersArray,
    *elementInnerDiametersArray,
    *elementBoundaryDiametersArray,
    *elementBarycentersArray,
    *elementBoundaryBarycentersArray,
    *nodeDiametersArray,
    *nodeSupportArray,
    *newestNodeBases,
    *elementIJK,             //NURBS
    *weights,                 //NURBS
    *U_KNOT, *V_KNOT,*W_KNOT; //NURBS
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  int dims[3]; 
  dims[0] = MESH(cmesh).nElements_global;
  dims[1] = MESH(cmesh).nNodes_element;
  elementNodesArray = PyArray_FromDimsAndData(2,
                                              dims, 
                                              PyArray_INT, 
                                              (char*)MESH(cmesh).elementNodesArray);
  if (MESH(cmesh).nodeElementOffsets)
    dims[0] = MESH(cmesh).nodeElementOffsets[MESH(cmesh).nNodes_global];
  else
    dims[0] = 0;
  nodeElementsArray = PyArray_FromDimsAndData(1,
                                              dims,
                                              PyArray_INT, 
                                              (char*)MESH(cmesh).nodeElementsArray);
  dims[0] = MESH(cmesh).nNodes_global+1;
  nodeElementOffsets = PyArray_FromDimsAndData(1,
					       dims,
					       PyArray_INT,
					       (char*)MESH(cmesh).nodeElementOffsets);
  dims[0] = MESH(cmesh).nElements_global;
  dims[1] = MESH(cmesh).nElementBoundaries_element;
  elementNeighborsArray = PyArray_FromDimsAndData(2,
                                                  dims,
                                                  PyArray_INT,
                                                  (char*)MESH(cmesh).elementNeighborsArray);
  dims[0] = MESH(cmesh).nElements_global;
  dims[1] = MESH(cmesh).nElementBoundaries_element;
  elementBoundariesArray = PyArray_FromDimsAndData(2,
						   dims,
						   PyArray_INT,
						   (char*)MESH(cmesh).elementBoundariesArray);
  dims[0] = MESH(cmesh).nElementBoundaries_global;
  dims[1] = MESH(cmesh).nNodes_elementBoundary;
  elementBoundaryNodesArray = PyArray_FromDimsAndData(2,
                                                      dims,
                                                      PyArray_INT,
                                                      (char*)MESH(cmesh).elementBoundaryNodesArray);
  dims[0] = MESH(cmesh).nElementBoundaries_global;
  dims[1] = 2;
  elementBoundaryElementsArray = PyArray_FromDimsAndData(2,
                                                         dims,
                                                         PyArray_INT,
                                                         (char*)MESH(cmesh).elementBoundaryElementsArray);
  dims[0] = MESH(cmesh).nElementBoundaries_global;
  dims[1] = 2;
  elementBoundaryLocalElementBoundariesArray = PyArray_FromDimsAndData(2,
                                                                       dims,
                                                                       PyArray_INT,
                                                                       (char*)MESH(cmesh).elementBoundaryLocalElementBoundariesArray);
  dims[0] = MESH(cmesh).nInteriorElementBoundaries_global;
  interiorElementBoundariesArray = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).interiorElementBoundariesArray);
  dims[0] = MESH(cmesh).nExteriorElementBoundaries_global;
  exteriorElementBoundariesArray = PyArray_FromDimsAndData(1,
                                                           dims,
                                                           PyArray_INT,
                                                           (char*)MESH(cmesh).exteriorElementBoundariesArray);
  dims[0] = MESH(cmesh).nEdges_global;
  dims[1] = 2;
  edgeNodesArray = PyArray_FromDimsAndData(2,
                                           dims,
                                           PyArray_INT,
                                           (char*)MESH(cmesh).edgeNodesArray);
  if (MESH(cmesh).nodeStarOffsets != NULL)
    dims[0] = MESH(cmesh).nodeStarOffsets[MESH(cmesh).nNodes_global];
  else
    dims[0] = 0;
  nodeStarArray = PyArray_FromDimsAndData(1,
                                          dims,
                                          PyArray_INT,
                                          (char*)MESH(cmesh).nodeStarArray);
  dims[0] = MESH(cmesh).nNodes_global+1;
  nodeStarOffsets = PyArray_FromDimsAndData(1,
                                            dims,
                                            PyArray_INT,
                                            (char*)MESH(cmesh).nodeStarOffsets);
  dims[0] = MESH(cmesh).nElements_global;
  elementMaterialTypes = PyArray_FromDimsAndData(1,
                                                dims,
                                                PyArray_INT,
                                                (char*)MESH(cmesh).elementMaterialTypes);
  dims[0] = MESH(cmesh).nElementBoundaries_global;
  elementBoundaryMaterialTypes = PyArray_FromDimsAndData(1,
							 dims,
							 PyArray_INT,
							 (char*)MESH(cmesh).elementBoundaryMaterialTypes);
  dims[0] = MESH(cmesh).nNodes_global;
  nodeMaterialTypes = PyArray_FromDimsAndData(1,
					      dims,
					      PyArray_INT,
					      (char*)MESH(cmesh).nodeMaterialTypes);
  dims[0] = MESH(cmesh).nNodes_global;
  dims[1] = 3;
  nodeArray = PyArray_FromDimsAndData(2,
                                     dims,
                                     PyArray_DOUBLE,
                                     (char*)MESH(cmesh).nodeArray);
  dims[0] = MESH(cmesh).nElements_global;
  elementDiametersArray = PyArray_FromDimsAndData(1,
                                      dims,
                                      PyArray_DOUBLE,
                                      (char*)MESH(cmesh).elementDiametersArray);
  dims[0] = MESH(cmesh).nElements_global;
  elementInnerDiametersArray = PyArray_FromDimsAndData(1,
                                      dims,
                                      PyArray_DOUBLE,
                                      (char*)MESH(cmesh).elementInnerDiametersArray);
  dims[0] = MESH(cmesh).nElementBoundaries_global;
  elementBoundaryDiametersArray = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_DOUBLE,
                                                          (char*)MESH(cmesh).elementBoundaryDiametersArray);
  dims[0] = MESH(cmesh).nElements_global;
  dims[1] = 3;
  elementBarycentersArray = PyArray_FromDimsAndData(2,
						    dims,
						    PyArray_DOUBLE,
						    (char*)MESH(cmesh).elementBarycentersArray);
  dims[0] = MESH(cmesh).nElementBoundaries_global;
  dims[1] = 3;
  elementBoundaryBarycentersArray = PyArray_FromDimsAndData(2,
							    dims,
							    PyArray_DOUBLE,
							    (char*)MESH(cmesh).elementBoundaryBarycentersArray);

  dims[0] = MESH(cmesh).nNodes_global;
  nodeDiametersArray = PyArray_FromDimsAndData(1,
					       dims,
					       PyArray_DOUBLE,
					       (char*)MESH(cmesh).nodeDiametersArray);
  
  dims[0] = MESH(cmesh).nNodes_global;
  nodeSupportArray = PyArray_FromDimsAndData(1,
					     dims,
					     PyArray_DOUBLE,
					     (char*)MESH(cmesh).nodeSupportArray);
  
// NURBS
  //mwf elementIJK,weights  V_KNOT,W_KNOT,U_KNOT unitialized? switched to MESH(cmesh).
  dims[0] = 0;
  if (MESH(cmesh).elementIJK != NULL) dims[0] = MESH(cmesh).nElements_global*3;
  elementIJK = PyArray_FromDimsAndData(1,
							    dims,
							    PyArray_DOUBLE,
							    (char*)MESH(cmesh).elementIJK);
  dims[0] = 0;
  if (MESH(cmesh).weights != NULL) dims[0] = MESH(cmesh).nElements_global;
  weights = PyArray_FromDimsAndData(1,
							    dims,
							    PyArray_DOUBLE,
							    (char*)MESH(cmesh).weights);
  dims[0] = 0;
  if (MESH(cmesh).U_KNOT != NULL) dims[0] = MESH(cmesh).nx+MESH(cmesh).px+1;
  U_KNOT = PyArray_FromDimsAndData(1,
							    dims,
							    PyArray_DOUBLE,
							    (char*)MESH(cmesh).U_KNOT);
  dims[0] = 0;
  if (MESH(cmesh).V_KNOT != NULL) dims[0] = MESH(cmesh).ny+MESH(cmesh).py+1;
  V_KNOT = PyArray_FromDimsAndData(1,
							    dims,
							    PyArray_DOUBLE,
							    (char*)MESH(cmesh).V_KNOT);
  dims[0] = 0;
  if (MESH(cmesh).W_KNOT != NULL) dims[0] = MESH(cmesh).nz+MESH(cmesh).pz+1;
  W_KNOT = PyArray_FromDimsAndData(1,
							    dims,
							    PyArray_DOUBLE,
							    (char*)MESH(cmesh).W_KNOT);							    
// NURBS

  return Py_BuildValue("iiiiiiiiiiiOOOOOOOOOOOOOOOOOiiiiiiOOOOOOOOOOOOdddd",
                       MESH(cmesh).nElements_global,
                       MESH(cmesh).nNodes_global,
                       MESH(cmesh).nNodes_element,
                       MESH(cmesh).nNodes_elementBoundary,
                       MESH(cmesh).nElementBoundaries_element,
                       MESH(cmesh).nElementBoundaries_global,
                       MESH(cmesh).nInteriorElementBoundaries_global,
                       MESH(cmesh).nExteriorElementBoundaries_global,
                       MESH(cmesh).max_nElements_node,
                       MESH(cmesh).nEdges_global,
                       MESH(cmesh).max_nNodeNeighbors_node,
                       elementNodesArray, //the nodes numbers for each element
                       nodeElementsArray,        //the element numbers for each node
                       nodeElementOffsets,       //offsets for indexing into nodeElementsArray
                       elementNeighborsArray,    //the elment numbers for each neighboring element
                       elementBoundariesArray,   //element boundaries for each element
                       elementBoundaryNodesArray,
                       elementBoundaryElementsArray, //the element numbers for each element boundary
                       elementBoundaryLocalElementBoundariesArray, //the local element boundary number for the left and right neighbors of the element boundary
                       interiorElementBoundariesArray, //the element boundary numbers of the interior element boundaries
                       exteriorElementBoundariesArray, //the elemetn boundary numbers of the exterior element boundaries
                       edgeNodesArray,
                       nodeStarArray,
                       nodeStarOffsets,
                       elementMaterialTypes,
                       elementBoundaryMaterialTypes,
                       nodeMaterialTypes,
                       nodeArray,
                       MESH(cmesh).nx,MESH(cmesh).ny,MESH(cmesh).nz,        //NURBS 
                       MESH(cmesh).px,MESH(cmesh).py,MESH(cmesh).pz,        //NURBS                       
                       elementIJK, //NURBS                 
                       weights,    //NURBS                   
                       U_KNOT,     //NURBS    
                       V_KNOT,     //NURBS    
                       W_KNOT,     //NURBS   
                       elementDiametersArray,
                       elementInnerDiametersArray,
                       elementBoundaryDiametersArray,
		       elementBarycentersArray,
                       elementBoundaryBarycentersArray,
		       nodeDiametersArray,
		       nodeSupportArray,
                       MESH(cmesh).h,
                       MESH(cmesh).hMin,
                       MESH(cmesh).sigmaMax,
                       MESH(cmesh).volume);
}

static PyObject* cmeshToolsBuildPythonMeshInterfaceNoArrays(PyObject* self,
							    PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  return Py_BuildValue("iiiiiiiiiiidddd",
                       MESH(cmesh).nElements_global,
                       MESH(cmesh).nNodes_global,
                       MESH(cmesh).nNodes_element,
                       MESH(cmesh).nNodes_elementBoundary,
                       MESH(cmesh).nElementBoundaries_element,
                       MESH(cmesh).nElementBoundaries_global,
                       MESH(cmesh).nInteriorElementBoundaries_global,
                       MESH(cmesh).nExteriorElementBoundaries_global,
                       MESH(cmesh).max_nElements_node,
                       MESH(cmesh).nEdges_global,
                       MESH(cmesh).max_nNodeNeighbors_node,
                       MESH(cmesh).h,
                       MESH(cmesh).hMin,
                       MESH(cmesh).sigmaMax,
                       MESH(cmesh).volume);
}
static PyObject* cmeshToolsBuildLevel0PythonMeshInterface(PyObject* self,
							  PyObject* args)
{
  PyObject *cmesh,
    *elementNodesArray,
    *elementMaterialTypes,
    *nodeMaterialTypes,
    *nodeArray;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  int dims[3]; 
  dims[0] = MESH(cmesh).nElements_global;
  dims[1] = MESH(cmesh).nNodes_element;
  elementNodesArray = PyArray_FromDimsAndData(2,
                                              dims, 
                                              PyArray_INT, 
                                              (char*)MESH(cmesh).elementNodesArray);
  dims[0] = MESH(cmesh).nElements_global;
  elementMaterialTypes = PyArray_FromDimsAndData(1,
                                                dims,
                                                PyArray_INT,
                                                (char*)MESH(cmesh).elementMaterialTypes);
  dims[0] = MESH(cmesh).nNodes_global;
  nodeMaterialTypes = PyArray_FromDimsAndData(1,
					      dims,
					      PyArray_INT,
					      (char*)MESH(cmesh).nodeMaterialTypes);
  dims[0] = MESH(cmesh).nNodes_global;
  dims[1] = 3;
  nodeArray = PyArray_FromDimsAndData(2,
                                     dims,
                                     PyArray_DOUBLE,
                                     (char*)MESH(cmesh).nodeArray);
  return Py_BuildValue("iiiOOOO",
                       MESH(cmesh).nElements_global,
                       MESH(cmesh).nNodes_global,
                       MESH(cmesh).nNodes_element,
                       elementNodesArray, //the nodes numbers for each element
                       elementMaterialTypes,
                       nodeMaterialTypes,
                       nodeArray);
}

typedef struct
{
  PyObject_HEAD
  MultilevelMesh multilevelMesh;
} CMultilevelMesh;

#define MULTILEVELMESH(p) ((CMultilevelMesh*)p)->multilevelMesh

static PyObject*
CMultilevelMesh_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  CMultilevelMesh *self;
  self = (CMultilevelMesh *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
CMultilevelMesh_init(CMultilevelMesh *self, PyObject *args, PyObject *kwds)
{
  int nLevels;
  PyObject *cmesh;
  if(!PyArg_ParseTuple(args,
                       "Oi",
                       &cmesh,
                       &nLevels))
    return -1;
  if (MESH(cmesh).nNodes_element == 2)
    {
      globallyRefineEdgeMesh(nLevels,MESH(cmesh),self->multilevelMesh);
      for (int i=1;i<nLevels;i++)
        {
          constructElementBoundaryElementsArray_edge(self->multilevelMesh.meshArray[i]);
          allocateGeometricInfo_edge(self->multilevelMesh.meshArray[i]);
          computeGeometricInfo_edge(self->multilevelMesh.meshArray[i]);
	  assignElementBoundaryMaterialTypesFromParent(self->multilevelMesh.meshArray[i-1],
						       self->multilevelMesh.meshArray[i],
						       self->multilevelMesh.elementParentsArray[i],
						       1);
        }
    }
  else if (MESH(cmesh).nNodes_element == 3)
    {
      globallyRefineTriangularMesh(nLevels,MESH(cmesh),self->multilevelMesh);
      for (int i=1;i<nLevels;i++)
        {
          constructElementBoundaryElementsArray_triangle(self->multilevelMesh.meshArray[i]);
          allocateGeometricInfo_triangle(self->multilevelMesh.meshArray[i]);
          computeGeometricInfo_triangle(self->multilevelMesh.meshArray[i]);
	  assignElementBoundaryMaterialTypesFromParent(self->multilevelMesh.meshArray[i-1],
						       self->multilevelMesh.meshArray[i],
						       self->multilevelMesh.elementParentsArray[i],
						       2);
        }
    }
  else
    {
      globallyRefineTetrahedralMesh(nLevels,MESH(cmesh),self->multilevelMesh);
      for (int i=1;i<nLevels;i++)
        {
          constructElementBoundaryElementsArray_tetrahedron(self->multilevelMesh.meshArray[i]);
          allocateGeometricInfo_tetrahedron(self->multilevelMesh.meshArray[i]);
          computeGeometricInfo_tetrahedron(self->multilevelMesh.meshArray[i]);
	  assignElementBoundaryMaterialTypesFromParent(self->multilevelMesh.meshArray[i-1],
						       self->multilevelMesh.meshArray[i],
						       self->multilevelMesh.elementParentsArray[i],
						       3);
        }
    }
  return 0;
}

static  void
CMultilevelMesh_dealloc(CMultilevelMesh* self)
{
  deleteMultilevelMesh(self->multilevelMesh);
}

static PyTypeObject CMultilevelMeshType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "cmeshTools.CMultilevelMesh",             /*tp_name*/
  sizeof(CMultilevelMesh), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)CMultilevelMesh_dealloc,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,        /*tp_flags*/
  "CMultilevelMesh objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  0,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)CMultilevelMesh_init,      /* tp_init */
  0,                         /* tp_alloc */
  CMultilevelMesh_new,                 /* tp_new */
};

static PyObject* cmeshToolsBuildPythonMultilevelMeshInterface(PyObject* self,
                                                              PyObject* args)
{
  PyObject *cmultilevelMesh,
    *mesh,
    *meshList,
    *elementChildrenArrayList,
    *elementChildrenOffsetsList,
    *elementParentsArrayList;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmultilevelMesh))
    return NULL;
  meshList = PyList_New(0);
  mesh = CMesh_FromMesh(&MULTILEVELMESH(cmultilevelMesh).meshArray[0]);
  PyList_Append(meshList,mesh);
  elementChildrenArrayList = PyList_New(0);
  elementChildrenOffsetsList = PyList_New(0);
  elementParentsArrayList  = PyList_New(0); //put in one empty entry
  int n,dims[1];
  npy_intp newdims[1];
  newdims[0]=0;
  PyList_Append(elementParentsArrayList,PyArray_SimpleNew(1,newdims,PyArray_INT));
  for(n=1;n<MULTILEVELMESH(cmultilevelMesh).nLevels;n++)
    {
      mesh = CMesh_FromMesh(&MULTILEVELMESH(cmultilevelMesh).meshArray[n]);
      PyList_Append(meshList,mesh);

      dims[0] = MULTILEVELMESH(cmultilevelMesh).meshArray[n].nElements_global;
      PyList_Append(elementParentsArrayList,
                    PyArray_FromDimsAndData(1,
                                            dims, 
                                            PyArray_INT, 
                                            (char*)(MULTILEVELMESH(cmultilevelMesh).elementParentsArray[n])));
      
      dims[0] = MULTILEVELMESH(cmultilevelMesh).elementChildrenOffsets[n-1][MULTILEVELMESH(cmultilevelMesh).meshArray[n-1].nElements_global];
      PyList_Append(elementChildrenArrayList,
                    PyArray_FromDimsAndData(1,
                                            dims, 
                                            PyArray_INT, 
                                            (char*)(MULTILEVELMESH(cmultilevelMesh).elementChildrenArray[n-1])));
      
      dims[0] = MULTILEVELMESH(cmultilevelMesh).meshArray[n-1].nElements_global+1;
      PyList_Append(elementChildrenOffsetsList,
                    PyArray_FromDimsAndData(1,
                                            dims, 
                                            PyArray_INT, 
                                            (char*)(MULTILEVELMESH(cmultilevelMesh).elementChildrenOffsets[n-1])));
    }
  return Py_BuildValue("iOOOO",
                       MULTILEVELMESH(cmultilevelMesh).nLevels,
                       meshList,
                       elementParentsArrayList,
                       elementChildrenArrayList,
                       elementChildrenOffsetsList);
}


static PyObject* cmeshToolsLocallyRefineMultilevelMesh(PyObject* self,
						       PyObject* args)
{
  PyObject *cmultilevelMesh,
    *elementTagArray;
  int nSpace,failed,finestLevel,refineTypeFlag;
  refineTypeFlag = 0;
  if (!PyArg_ParseTuple(args,
			"iOO|i",			
			&nSpace,
			&cmultilevelMesh,
			&elementTagArray,
			&refineTypeFlag))
    return NULL;

  if (nSpace == 1)
    {
      failed = locallyRefineEdgeMesh(MULTILEVELMESH(cmultilevelMesh),
				     IDATA(elementTagArray));
      finestLevel = MULTILEVELMESH(cmultilevelMesh).nLevels;
      constructElementBoundaryElementsArray_edge(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1]);
      allocateGeometricInfo_edge(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1]);
      computeGeometricInfo_edge(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1]);
      if (finestLevel > 1)
	assignElementBoundaryMaterialTypesFromParent(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-2],
						     MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1],
						     MULTILEVELMESH(cmultilevelMesh).elementParentsArray[finestLevel-1],
						     1);
    }
  else if (nSpace == 2)
    {
      if (refineTypeFlag==1)
	{
	  failed = locallyRefineTriangleMesh_4T(MULTILEVELMESH(cmultilevelMesh),
						IDATA(elementTagArray));
	}
      else if (refineTypeFlag==2)
        {
          failed = locallyRefineTriangleMesh_redGreen(MULTILEVELMESH(cmultilevelMesh),IDATA(elementTagArray));
        }
      else
	{
	  failed = locallyRefineTriangleMesh(MULTILEVELMESH(cmultilevelMesh),
					     IDATA(elementTagArray));
	}
      finestLevel = MULTILEVELMESH(cmultilevelMesh).nLevels;
      //std::cout<<"eb info"<<std::endl;
      constructElementBoundaryElementsArray_triangle(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1]);
      //std::cout<<"geo info"<<std::endl;
      allocateGeometricInfo_triangle(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1]);
      computeGeometricInfo_triangle(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1]);
      //std::cout<<"eb types"<<std::endl;
      if (finestLevel > 1)
	assignElementBoundaryMaterialTypesFromParent(MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-2],
						     MULTILEVELMESH(cmultilevelMesh).meshArray[finestLevel-1],
						     MULTILEVELMESH(cmultilevelMesh).elementParentsArray[finestLevel-1],
						     2);

    }
  else
    {
      std::cout<<"locallyRefine nSpace= "<<nSpace<<" not implemented! Returning"<<std::endl;
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cmeshToolsSetNewestNodeBasesToLongestEdge(PyObject* self,
							   PyObject* args)
{
  //PyObject *cmesh;
  PyObject *cmultilevelMesh;
  int nSpace,failed;
  if (!PyArg_ParseTuple(args,
			"iO",			
			&nSpace,
			&cmultilevelMesh))
    return NULL;

  if (nSpace == 2)
    {
      failed = setNewestNodeBasesToLongestEdge(MULTILEVELMESH(cmultilevelMesh));
    }
  else
    {
      std::cout<<"setNewestNodeBases= "<<nSpace<<" not implemented! Returning"<<std::endl;
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cmeshToolsGenerateTetrahedralMeshFromRectangularGrid(PyObject* self,
                                                                      PyObject* args)
{
  int nx,ny,nz;
  double Lx,Ly,Lz;
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "iiidddO",
                        &nx,
                        &ny,
                        &nz,
                        &Lx,
                        &Ly,
                        &Lz,
                        &cmesh))
    return NULL;
  regularHexahedralToTetrahedralMeshElements(nx,ny,nz,MESH(cmesh));
  regularHexahedralToTetrahedralMeshNodes(nx,ny,nz,Lx,Ly,Lz,MESH(cmesh));
  constructElementBoundaryElementsArray_tetrahedron(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsComputeGeometricInfo_tetrahedron(PyObject* self,
                                                            PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  computeGeometricInfo_tetrahedron(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsAllocateGeometricInfo_tetrahedron(PyObject* self,
							   PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  allocateGeometricInfo_tetrahedron(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsAllocateNodeAndElementNodeDataStructures(PyObject* self,
								    PyObject* args)
{
  PyObject *cmesh;
  int nElements_global, nNodes_global, nNodes_element;
  if (!PyArg_ParseTuple(args,
                        "Oiii",
                        &cmesh,
			&nElements_global,
			&nNodes_global,
			&nNodes_element))
    return NULL;
  allocateNodeAndElementNodeDataStructures(MESH(cmesh),nElements_global,nNodes_global,nNodes_element);
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cmeshToolsConstructElementBoundaryElementsArray(PyObject* self,
								 PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  if (MESH(cmesh).nNodes_element == 4)
    {
      constructElementBoundaryElementsArray_tetrahedron(MESH(cmesh));
    }
  else if (MESH(cmesh).nNodes_element == 3)
    {
      constructElementBoundaryElementsArray_triangle(MESH(cmesh));
    }
  else
    {
      constructElementBoundaryElementsArray_edge(MESH(cmesh));
    }
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsGenerateTriangularMeshFromRectangularGrid(PyObject* self,
                                                                     PyObject* args)
{
  int nx,ny;
  double Lx,Ly;
  PyObject *cmesh;
  int triangleFlag=1;
  if (!PyArg_ParseTuple(args,
                        "iiddO|i",
                        &nx,
                        &ny,
                        &Lx,
                        &Ly,
                        &cmesh,
			&triangleFlag))
    return NULL;
  regularRectangularToTriangularMeshElements(nx,ny,MESH(cmesh),triangleFlag);
  regularRectangularToTriangularMeshNodes(nx,ny,Lx,Ly,MESH(cmesh));
  constructElementBoundaryElementsArray_triangle(MESH(cmesh));
  regularRectangularToTriangularElementBoundaryMaterials(Lx,Ly,MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsGenerateHexahedralMeshFromRectangularGrid(PyObject* self,
                                                                     PyObject* args)
{
  int nx,ny,nz,px,py,pz;
  double Lx,Ly,Lz;
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "iiiiiidddO",
                        &nx,
                        &ny,
                        &nz,
			&px,
			&py,
			&pz,
                        &Lx,
                        &Ly,
                        &Lz,
                        &cmesh))
    return NULL;
  regularHexahedralMeshElements(nx,ny,nz,px,py,pz,MESH(cmesh));
  regularMeshNodes(nx,ny,nz,Lx,Ly,Lz,MESH(cmesh));
  constructElementBoundaryElementsArray_hexahedron(MESH(cmesh));
  regularHexahedralToTetrahedralElementBoundaryMaterials(Lx,Ly,Lz,MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

#define TRIANGULATEIO(p) ((triangulateio*) p)

static PyObject* 
cmeshToolsGenerateFromTriangleMesh(PyObject* self,
				    PyObject* args)
{
  PyObject *cmesh,*cob;
  triangulateio* trimesh;
  int base;
  if (!PyArg_ParseTuple(args,
                        "OOi",
                        &cmesh,
                        &cob,
			&base))
    return NULL;
  trimesh = (triangulateio *) PyCObject_AsVoidPtr(cob);

  setFromTriangleElements(TRIANGULATEIO(trimesh),MESH(cmesh),base);
  setFromTriangleNodes(TRIANGULATEIO(trimesh),MESH(cmesh),base);
  constructElementBoundaryElementsArray_triangle(MESH(cmesh));
  copyElementBoundaryMaterialTypesFromTriangle(TRIANGULATEIO(trimesh),
					       MESH(cmesh),base);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* 
cmeshToolsGenerateFromTriangleFiles(PyObject* self,
				    PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = readTriangleMesh(MESH(cmesh),filebase,base);
  constructElementBoundaryElementsArray_triangle(MESH(cmesh));
  failed = readTriangleElementBoundaryMaterialTypes(MESH(cmesh),filebase,base);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* 
cmeshToolsWriteTriangleFiles(PyObject* self,
			     PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = writeTriangleMesh(MESH(cmesh),filebase,base);

  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* 
cmeshToolsGenerateFromTetgenFiles(PyObject* self,
				  PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = readTetgenMesh(MESH(cmesh),filebase,base);
  constructElementBoundaryElementsArray_tetrahedron(MESH(cmesh));
  failed = readTetgenElementBoundaryMaterialTypes(MESH(cmesh),filebase,base);

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* 
cmeshToolsGenerateFromTetgenFilesParallel(PyObject* self,
					  PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = readTetgenMesh(MESH(cmesh),filebase,base);
  constructElementBoundaryElementsArray_tetrahedron(MESH(cmesh));
  failed = readTetgenElementBoundaryMaterialTypes(MESH(cmesh),filebase,base);

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* 
cmeshToolsWriteTetgenFiles(PyObject* self,
			   PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = writeTetgenMesh(MESH(cmesh),filebase,base);

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* 
cmeshToolsWrite3dmFiles(PyObject* self,
			   PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = write3dmMesh(MESH(cmesh),filebase,base);

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* 
cmeshToolsWrite2dmFiles(PyObject* self,
			   PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = write2dmMesh(MESH(cmesh),filebase,base);

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* 
cmeshToolsGenerateFromHexFile(PyObject* self,
			      PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = readHex(MESH(cmesh),filebase,base);
  constructElementBoundaryElementsArray_hexahedron(MESH(cmesh));
  //failed = readBC(MESH(cmesh),filebase,base);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* 
cmeshToolsGenerateFrom3DMFile(PyObject* self,
			      PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = read3DM(MESH(cmesh),filebase,base);
  constructElementBoundaryElementsArray_tetrahedron(MESH(cmesh));
  failed = readBC(MESH(cmesh),filebase,base);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* 
cmeshToolsGenerateFrom2DMFile(PyObject* self,
			      PyObject* args)
{
  PyObject *cmesh;
  const char *filebase;
  int base,failed;
  if (!PyArg_ParseTuple(args,
                        "Osi",
                        &cmesh,
                        &filebase,
			&base))
    return NULL;

  failed = read2DM(MESH(cmesh),filebase,base);
  constructElementBoundaryElementsArray_triangle(MESH(cmesh));
  failed = readBC(MESH(cmesh),filebase,base);
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* cmeshToolsComputeGeometricInfo_triangle(PyObject* self,
                                                         PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  computeGeometricInfo_triangle(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsAllocateGeometricInfo_triangle(PyObject* self,
							PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  allocateGeometricInfo_triangle(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsGenerateEdgeMeshFromRectangularGrid(PyObject* self,
                                                               PyObject* args)
{
  int nx;
  double Lx;
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "idO",
                        &nx,
                        &Lx,
                        &cmesh))
    return NULL;
  edgeMeshElements(nx,MESH(cmesh));
  regularEdgeMeshNodes(nx,Lx,MESH(cmesh));
  constructElementBoundaryElementsArray_edge(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsComputeGeometricInfo_edge(PyObject* self,
                                                     PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  computeGeometricInfo_edge(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsAllocateGeometricInfo_edge(PyObject* self,
                                                     PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  allocateGeometricInfo_edge(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject* cmeshToolsComputeGeometricInfo_hexahedron(PyObject* self,
                                                     PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  computeGeometricInfo_hexahedron(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsAllocateGeometricInfo_hexahedron(PyObject* self,
                                                     PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  allocateGeometricInfo_hexahedron(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsComputeGeometricInfo_NURBS(PyObject* self,
                                                     PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  computeGeometricInfo_NURBS(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsAllocateGeometricInfo_NURBS(PyObject* self,
                                                      PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  allocateGeometricInfo_NURBS(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}


typedef struct
{
  PyObject_HEAD
  std::map<int,std::set<int> > columnIndecesMap;
  std::map<int,std::map<int,int> > columnOffsetsMap;
} SparsityInfo;

extern "C"
{
static PyObject*
SparsityInfo_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  SparsityInfo *self;
  self = (SparsityInfo *)type->tp_alloc(type,0);
  return (PyObject*)self;
}

static int
SparsityInfo_init(SparsityInfo *self, PyObject *args, PyObject *kwds)
{
  new(&self->columnIndecesMap) std::map<int,std::set<int> >;
  new(&self->columnOffsetsMap) std::map<int,std::map<int,int> >;
  return 0;
}

static void
SparsityInfo_dealloc(SparsityInfo *self)
{
  self->columnIndecesMap.clear();
  self->columnOffsetsMap.clear();
  self->ob_type->tp_free((PyObject*)self);
}


static PyObject* SparsityInfo_findNonzeros(SparsityInfo *self,
                                           PyObject *args)
{
  int nElements_global,nDOF_test_element,nDOF_trial_element;
  PyObject *Py_nFreeDOF_test,*Py_freeGlobal_test,*Py_nFreeDOF_trial,*Py_freeGlobal_trial;
  int offset_test,stride_test,offset_trial,stride_trial;
  int hasNumericalFlux,hasDiffusionInMixedForm,needNumericalFluxJacobian;
  int nElementBoundaries_element;
  PyObject *Py_elementNeighborsArray;
  int nInteriorElementBoundaries_global;
  PyObject *Py_interiorElementBoundariesArray;
  PyObject *Py_elementBoundaryElementsArray;
  PyObject *Py_elementBoundaryLocalElementBoundariesArray;
  int hasFluxBoundaryConditions;
  int nExteriorElementBoundaries_global;
  PyObject *Py_exteriorElementBoundariesArray;
  int hasOutflowBoundary,needOutflowJacobian;
  if(!PyArg_ParseTuple(args,
                       "iiiOOOOiiiiiiiiOiOOOiiOii",
                       &nElements_global,
                       &nDOF_test_element,
                       &nDOF_trial_element,
                       &Py_nFreeDOF_test,
                       &Py_freeGlobal_test,
                       &Py_nFreeDOF_trial,
                       &Py_freeGlobal_trial,
                       &offset_test,
                       &stride_test,
                       &offset_trial,
                       &stride_trial,
                       &hasNumericalFlux,
                       &hasDiffusionInMixedForm,
                       &needNumericalFluxJacobian,
                       &nElementBoundaries_element,
                       &Py_elementNeighborsArray,
                       &nInteriorElementBoundaries_global,
                       &Py_interiorElementBoundariesArray,
                       &Py_elementBoundaryElementsArray,
                       &Py_elementBoundaryLocalElementBoundariesArray,
                       &hasFluxBoundaryConditions,
                       &nExteriorElementBoundaries_global,
                       &Py_exteriorElementBoundariesArray,
                       &hasOutflowBoundary,
                       &needOutflowJacobian))
    return NULL;
  int *nFreeDOF_test=IDATA(Py_nFreeDOF_test),
    *freeGlobal_test=IDATA(Py_freeGlobal_test),
    *nFreeDOF_trial=IDATA(Py_nFreeDOF_trial),
    *freeGlobal_trial=IDATA(Py_freeGlobal_trial),
    *elementNeighborsArray=IDATA(Py_elementNeighborsArray),
    *interiorElementBoundariesArray=IDATA(Py_interiorElementBoundariesArray),
    *elementBoundaryElementsArray=IDATA(Py_elementBoundaryElementsArray),
    *elementBoundaryLocalElementBoundariesArray=IDATA(Py_elementBoundaryLocalElementBoundariesArray),
    *exteriorElementBoundariesArray=IDATA(Py_exteriorElementBoundariesArray);
  //elements
  for(int eN=0;eN<nElements_global;eN++)
    {
      for(int ii=0;ii<nFreeDOF_test[eN];ii++)
        {
          int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
          for(int jj=0;jj<nFreeDOF_trial[eN];jj++)
            {
              int J = offset_trial + stride_trial*freeGlobal_trial[eN*nDOF_trial_element+jj];
              self->columnIndecesMap[I].insert(J);
            }
        }
      //get element neighbor DOF for mixed form diffusion
      if (hasNumericalFlux &&
          hasDiffusionInMixedForm &&
          needNumericalFluxJacobian)
        {
          for(int ebN=0;ebN<nElementBoundaries_element;ebN++)
            {
              int eN_ebN = elementNeighborsArray[eN*nElementBoundaries_element + ebN];
              if (eN_ebN >= 0)
                {
                  for (int ii=0;ii<nFreeDOF_test[eN];ii++)
                    {
                      int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
                      for (int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                        {
                          int J = offset_trial + stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                          self->columnIndecesMap[I].insert(J);
                        }
                    }
                }
            }
        }
    }
  //interior element boundaries
  if (hasNumericalFlux &&
      needNumericalFluxJacobian)
    {
      for(int ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
        {
          int ebN = interiorElementBoundariesArray[ebNI],
            left_eN_global   = elementBoundaryElementsArray[ebN*2 + 0],
            right_eN_global  = elementBoundaryElementsArray[ebN*2 + 1];
          //left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2 + 0],
          //right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2 + 1];
          for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
            {
              int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element + ii];
              for (int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                {
                  int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element + jj];
                  self->columnIndecesMap[left_I].insert(left_J);
                }
              for (int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                {
                  int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element + jj];
                  self->columnIndecesMap[left_I].insert(right_J);
                }
            }
          for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
            {
              int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
              for(int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                {
                  int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element+jj];
                  self->columnIndecesMap[right_I].insert(left_J);
                }
              for(int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                {
                  int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element+jj];
                  self->columnIndecesMap[right_I].insert(right_J);
                }
            }
          if(hasDiffusionInMixedForm)
            {
              for(int ebN_eN=0;ebN_eN<nElementBoundaries_element;ebN_eN++)
                {
                  int left_eN_ebN = elementNeighborsArray[left_eN_global*nElementBoundaries_element+ebN_eN],
                    right_eN_ebN = elementNeighborsArray[right_eN_global*nElementBoundaries_element+ebN_eN];
                  for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
                    {
                      int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element+ii];
                      if(left_eN_ebN >= 0)
                        {
                          for (int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                            {
                              int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                              self->columnIndecesMap[left_I].insert(left_J);
                            }
                        }
                      if(right_eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                            {
                              int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                              self->columnIndecesMap[left_I].insert(right_J);
                            }
                        }
                    }
                  for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
                    {
                      int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
                      if(left_eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                            {
                              int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                              self->columnIndecesMap[right_I].insert(left_J);
                            }
                        }
                      if(right_eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                            {
                              int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                              self->columnIndecesMap[right_I].insert(right_J);
                            }
                        }
                    }
                }
            }
        }
    }
  //exterior element boundaries
  if ((hasNumericalFlux &&
       needNumericalFluxJacobian) ||
      (hasOutflowBoundary &&
       needOutflowJacobian))
    {
      for(int ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
        {
          int ebN = exteriorElementBoundariesArray[ebNE],
            eN_global   = elementBoundaryElementsArray[ebN*2+0];
          //ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
          for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
            {
              int I = offset_test+stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
              for (int jj=0;jj<nFreeDOF_trial[eN_global];jj++)
                {
                  int J = offset_trial+stride_trial*freeGlobal_trial[eN_global*nDOF_trial_element+jj];
                  self->columnIndecesMap[I].insert(J);
                }
            }
          if(hasNumericalFlux &&
             hasDiffusionInMixedForm)
            {
              for(int ebN_eN=0;ebN_eN < nElementBoundaries_element;ebN_eN++)
                {
                  int eN_ebN = elementNeighborsArray[eN_global*nElementBoundaries_element+ebN_eN];
                  for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
                    {
                      int I = offset_test + stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
                      if(eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                            {
                              int J = offset_trial+stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                              self->columnIndecesMap[I].insert(J);
                            }
                        }
                    }
                }
            }
        }
    }
  //debug
//   for (std::map<int, std::set<int> >::iterator mit = self->columnIndecesMap.begin();mit != self->columnIndecesMap.end();mit++)
//     {
//       for(std::set<int>::iterator sit = mit->second.begin();sit!=mit->second.end();sit++)
//         std::cout<<*sit<<'\t';
//       std::cout<<std::endl;
//     }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* SparsityInfo_getOffsets_CSR(SparsityInfo *self,
                                           PyObject *args)
{
  int nElements_global,nDOF_test_element,nDOF_trial_element;
  PyObject *Py_nFreeDOF_test,*Py_freeGlobal_test,*Py_nFreeDOF_trial,*Py_freeGlobal_trial;
  int offset_test,stride_test,offset_trial,stride_trial;
  int hasNumericalFlux,hasDiffusionInMixedForm,needNumericalFluxJacobian;
  int nElementBoundaries_element;
  PyObject *Py_elementNeighborsArray;
  int nInteriorElementBoundaries_global;
  PyObject *Py_interiorElementBoundariesArray;
  PyObject *Py_elementBoundaryElementsArray;
  PyObject *Py_elementBoundaryLocalElementBoundariesArray;
  int hasFluxBoundaryConditions;
  int nExteriorElementBoundaries_global;
  PyObject *Py_exteriorElementBoundariesArray;
  int hasOutflowBoundary,needOutflowJacobian;
  PyObject *Py_rowptr,
    *Py_csrRowIndeces,
    *Py_csrColumnOffsets,
    *Py_csrColumnOffsets_eNebN,
    *Py_csrColumnOffsets_eb,
    *Py_csrColumnOffsets_eb_eNebN;
  if(!PyArg_ParseTuple(args,
                       "iiiOOOOiiiiiiiiOiOOOiiOiiOOOOOO",
                       &nElements_global,
                       &nDOF_test_element,
                       &nDOF_trial_element,
                       &Py_nFreeDOF_test,
                       &Py_freeGlobal_test,
                       &Py_nFreeDOF_trial,
                       &Py_freeGlobal_trial,
                       &offset_test,
                       &stride_test,
                       &offset_trial,
                       &stride_trial,
                       &hasNumericalFlux,
                       &hasDiffusionInMixedForm,
                       &needNumericalFluxJacobian,
                       &nElementBoundaries_element,
                       &Py_elementNeighborsArray,
                       &nInteriorElementBoundaries_global,
                       &Py_interiorElementBoundariesArray,
                       &Py_elementBoundaryElementsArray,
                       &Py_elementBoundaryLocalElementBoundariesArray,
                       &hasFluxBoundaryConditions,
                       &nExteriorElementBoundaries_global,
                       &Py_exteriorElementBoundariesArray,
                       &hasOutflowBoundary,
                       &needOutflowJacobian,
                       &Py_rowptr,
                       &Py_csrRowIndeces,
                       &Py_csrColumnOffsets,
                       &Py_csrColumnOffsets_eNebN,
                       &Py_csrColumnOffsets_eb,
                       &Py_csrColumnOffsets_eb_eNebN))
    return NULL;
  int *nFreeDOF_test=IDATA(Py_nFreeDOF_test),
    *freeGlobal_test=IDATA(Py_freeGlobal_test),
    *nFreeDOF_trial=IDATA(Py_nFreeDOF_trial),
    *freeGlobal_trial=IDATA(Py_freeGlobal_trial),
    *elementNeighborsArray=IDATA(Py_elementNeighborsArray),
    *interiorElementBoundariesArray=IDATA(Py_interiorElementBoundariesArray),
    *elementBoundaryElementsArray=IDATA(Py_elementBoundaryElementsArray),
    *elementBoundaryLocalElementBoundariesArray=IDATA(Py_elementBoundaryLocalElementBoundariesArray),
    *exteriorElementBoundariesArray=IDATA(Py_exteriorElementBoundariesArray),
    *rowptr = IDATA(Py_rowptr),
    *csrRowIndeces=IDATA(Py_csrRowIndeces),
    *csrColumnOffsets=IDATA(Py_csrColumnOffsets),
    *csrColumnOffsets_eNebN=IDATA(Py_csrColumnOffsets_eNebN),
    *csrColumnOffsets_eb=IDATA(Py_csrColumnOffsets_eb),
    *csrColumnOffsets_eb_eNebN=IDATA(Py_csrColumnOffsets_eb_eNebN);
  //elements
  for(int eN=0;eN<nElements_global;eN++)
    {
      for(int ii=0;ii<nFreeDOF_test[eN];ii++)
        {
          int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
          csrRowIndeces[eN*nDOF_test_element+ii] = rowptr[I];
          for(int jj=0;jj<nFreeDOF_trial[eN];jj++)
            {
              int J = offset_trial + stride_trial*freeGlobal_trial[eN*nDOF_trial_element+jj];
              csrColumnOffsets[eN*nDOF_test_element*nDOF_trial_element+
                               ii*nDOF_trial_element+
                               jj] = self->columnOffsetsMap[I][J];
            }
        }
      //get element neighbor DOF for mixed form diffusion
      if (hasNumericalFlux &&
          hasDiffusionInMixedForm &&
          needNumericalFluxJacobian)
        {
          for(int ebN=0;ebN<nElementBoundaries_element;ebN++)
            {
              int eN_ebN = elementNeighborsArray[eN*nElementBoundaries_element + ebN];
              if (eN_ebN >= 0)
                {
                  for (int ii=0;ii<nFreeDOF_test[eN];ii++)
                    {
                      int I = offset_test + stride_test*freeGlobal_test[eN*nDOF_test_element+ii];
                      for (int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                        {
                          int J = offset_trial + stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                          csrColumnOffsets_eNebN[eN*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                 ebN*nDOF_test_element*nDOF_trial_element+
                                                 ii*nDOF_test_element+
                                                 jj] = self->columnOffsetsMap[I][J];
                        }
                    }
                }
            }
        }
    }
  //interior element boundaries
  if (hasNumericalFlux &&
      needNumericalFluxJacobian)
    {
      for(int ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
        {
          int ebN = interiorElementBoundariesArray[ebNI],
            left_eN_global   = elementBoundaryElementsArray[ebN*2 + 0],
            right_eN_global  = elementBoundaryElementsArray[ebN*2 + 1];
          //left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2 + 0],
          //right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2 + 1];
          for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
            {
              int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element + ii];
              for (int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                {
                  int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element + jj];
                  csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                      0*2*nDOF_test_element*nDOF_trial_element+
                                      0*nDOF_test_element*nDOF_trial_element+
                                      ii*nDOF_trial_element+
                                      jj] = self->columnOffsetsMap[left_I][left_J];
                }
              for (int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                {
                  int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element + jj];
                  csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                      0*2*nDOF_test_element*nDOF_trial_element+
                                      1*nDOF_test_element*nDOF_trial_element+
                                      ii*nDOF_trial_element+
                                      jj] = self->columnOffsetsMap[left_I][right_J];
                }
            }
          for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
            {
              int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
              for(int jj=0;jj<nFreeDOF_trial[left_eN_global];jj++)
                {
                  int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_global*nDOF_trial_element+jj];
                  csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                      1*2*nDOF_test_element*nDOF_trial_element+
                                      0*nDOF_test_element*nDOF_trial_element+
                                      ii*nDOF_trial_element+
                                      jj] = self->columnOffsetsMap[right_I][left_J];
                }
              for(int jj=0;jj<nFreeDOF_trial[right_eN_global];jj++)
                {
                  int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_global*nDOF_trial_element+jj];
                  csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                      1*2*nDOF_test_element*nDOF_trial_element+
                                      1*nDOF_test_element*nDOF_trial_element+
                                      ii*nDOF_trial_element+
                                      jj] = self->columnOffsetsMap[right_I][right_J];
                }
            }
          if(hasDiffusionInMixedForm)
            {
              for(int ebN_eN=0;ebN_eN<nElementBoundaries_element;ebN_eN++)
                {
                  int left_eN_ebN = elementNeighborsArray[left_eN_global*nElementBoundaries_element+ebN_eN],
                    right_eN_ebN = elementNeighborsArray[right_eN_global*nElementBoundaries_element+ebN_eN];
                  for(int ii=0;ii<nFreeDOF_test[left_eN_global];ii++)
                    {
                      int left_I = offset_test+stride_test*freeGlobal_test[left_eN_global*nDOF_test_element+ii];
                      if(left_eN_ebN >= 0)
                        {
                          for (int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                            {
                              int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                              csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        0*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        0*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                        ii*nDOF_trial_element+
                                                        jj] = self->columnOffsetsMap[left_I][left_J];
                            }
                        }
                      if(right_eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                            {
                              int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                              csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        0*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        1*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                        ii*nDOF_trial_element+
                                                        jj] = self->columnOffsetsMap[left_I][right_J];
                            }
                        }
                    }
                  for(int ii=0;ii<nFreeDOF_test[right_eN_global];ii++)
                    {
                      int right_I = offset_test+stride_test*freeGlobal_test[right_eN_global*nDOF_test_element+ii];
                      if(left_eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[left_eN_ebN];jj++)
                            {
                              int left_J = offset_trial+stride_trial*freeGlobal_trial[left_eN_ebN*nDOF_trial_element+jj];
                              csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        1*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        0*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                        ii*nDOF_trial_element+
                                                        jj] = self->columnOffsetsMap[right_I][left_J];
                            }
                        }
                      if(right_eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[right_eN_ebN];jj++)
                            {
                              int right_J = offset_trial+stride_trial*freeGlobal_trial[right_eN_ebN*nDOF_trial_element+jj];
                              csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        1*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        1*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                        ii*nDOF_trial_element+
                                                        jj] = self->columnOffsetsMap[right_I][right_J];
                            }
                        }
                    }
                }
            }
        }
    }
  //exterior element boundaries
  if ((hasNumericalFlux &&
       needNumericalFluxJacobian) ||
      (hasOutflowBoundary &&
       needOutflowJacobian))
    {
      for(int ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
        {
          int ebN = exteriorElementBoundariesArray[ebNE],
            eN_global   = elementBoundaryElementsArray[ebN*2+0];
          //ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
          for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
            {
              int I = offset_test+stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
              for (int jj=0;jj<nFreeDOF_trial[eN_global];jj++)
                {
                  int J = offset_trial+stride_trial*freeGlobal_trial[eN_global*nDOF_trial_element+jj];
                  csrColumnOffsets_eb[ebN*2*2*nDOF_test_element*nDOF_trial_element+
                                      0*2*nDOF_test_element*nDOF_trial_element+
                                      0*nDOF_test_element*nDOF_trial_element+
                                      ii*nDOF_trial_element+
                                      jj] = self->columnOffsetsMap[I][J];
                }
            }
          if(hasNumericalFlux &&
             hasDiffusionInMixedForm)
            {
              for(int ebN_eN=0;ebN_eN < nElementBoundaries_element;ebN_eN++)
                {
                  int eN_ebN = elementNeighborsArray[eN_global*nElementBoundaries_element+ebN_eN];
                  for (int ii=0;ii<nFreeDOF_test[eN_global];ii++)
                    {
                      int I = offset_test + stride_test*freeGlobal_test[eN_global*nDOF_test_element+ii];
                      if(eN_ebN >= 0)
                        {
                          for(int jj=0;jj<nFreeDOF_trial[eN_ebN];jj++)
                            {
                              int J = offset_trial+stride_trial*freeGlobal_trial[eN_ebN*nDOF_trial_element+jj];
                              csrColumnOffsets_eb_eNebN[ebN*2*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        0*2*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        0*nElementBoundaries_element*nDOF_test_element*nDOF_trial_element+
                                                        ebN_eN*nDOF_test_element*nDOF_trial_element+
                                                        ii*nDOF_trial_element+
                                                        jj] = self->columnOffsetsMap[I][J];
                            }
                        }
                    }
                }
            }
        }
    }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* SparsityInfo_getCSR(SparsityInfo *self,
                                     PyObject *args)
{
  //debug
//   for (std::map<int, std::set<int> >::iterator mit = self->columnIndecesMap.begin();mit != self->columnIndecesMap.end();mit++)
//     {
//       for(std::set<int>::iterator sit = mit->second.begin();sit!=mit->second.end();sit++)
//         std::cout<<*sit<<'\t';
//       std::cout<<std::endl;
//     }
  int nnz=0;
  int dim[1];
  dim[0] = int(self->columnIndecesMap.size())+1;
  PyArrayObject *rowptr = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_INT);
  int* rowptr_ = IDATA(rowptr);
  rowptr_[0] = 0;
  for(int I=1;I< self->columnIndecesMap.size()+1;I++)
    rowptr_[I]=rowptr_[I-1] + int(self->columnIndecesMap[I-1].size());
  nnz = rowptr_[self->columnIndecesMap.size()];
  dim[0] = nnz;
  PyArrayObject *colind = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_INT);
  PyArrayObject *nzval = (PyArrayObject *)PyArray_FromDims(1,dim,PyArray_DOUBLE);
  int *colind_ = IDATA(colind);
  int max_nonzeros=0;
  for(int I=0;I< self->columnIndecesMap.size();I++)
    {
      int offset=0;
      for(std::set<int>::iterator sit=self->columnIndecesMap[I].begin();sit != self->columnIndecesMap[I].end();sit++)
        {
          self->columnOffsetsMap[I][*sit] = offset;
          assert(rowptr_[I]+offset < nnz);
          colind_[rowptr_[I]+offset]=*sit;
          offset++;
        }
      std::sort(&colind_[rowptr_[I]],&colind_[rowptr_[I+1]]);
      max_nonzeros = std::max(max_nonzeros,rowptr_[I+1] - rowptr_[I]);
    }
  //std::cout<<"Proteus: Maximum nonzeros in any row is "<<max_nonzeros<<std::endl;
  return Py_BuildValue("(N,N,i,N)",PyArray_Return(rowptr),PyArray_Return(colind),nnz,PyArray_Return(nzval));
}


static PyObject* cmeshToolsGenerateNURBSMeshFromRectangularGrid(PyObject* self,
                                                                PyObject* args)
{
  int nx,ny,nz;
  int px,py,pz;
  double Lx,Ly,Lz;
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "iiiiiidddO",
                        &nx,
                        &ny,
                        &nz,
                        &px,
                        &py,
                        &pz,
                        &Lx,
                        &Ly,
                        &Lz,
                        &cmesh))
    return NULL;
  regularNURBSMeshElements(nx+px+1,ny+py+1,nz+pz+1,px,py,pz,MESH(cmesh));
  regularMeshNodes(nx+px+1,ny+py+1,nz+pz+1,Lx,Ly,Lz,MESH(cmesh));

  constructElementBoundaryElementsArray_NURBS(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cmeshToolsDeleteMeshDataStructures(PyObject* self,
						    PyObject* args)
{
  PyObject *cmesh;
  if (!PyArg_ParseTuple(args,
                        "O",
                        &cmesh))
    return NULL;
  deleteMesh(MESH(cmesh));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyMethodDef SparsityInfo_methods[] = {
  {"findNonzeros", 
   (PyCFunction)SparsityInfo_findNonzeros,
   METH_VARARGS, 
   "add nonzero entries to sparsity data structure"},
  {"getCSR", 
   (PyCFunction)SparsityInfo_getCSR,
   METH_NOARGS, 
   "get basic CSR data structures"},
  {"setOffsets_CSR", 
   (PyCFunction)SparsityInfo_getOffsets_CSR,
   METH_VARARGS, 
   "get offset arrys for indexing into CSR matrix"},
  {NULL}
};

static PyTypeObject SparsityInfoType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "cmeshTools.SparsityInfo",             /*tp_name*/
  sizeof(SparsityInfo), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor) SparsityInfo_dealloc,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,        /*tp_flags*/
  "SparsityInfo objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  SparsityInfo_methods,                         /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)SparsityInfo_init,      /* tp_init */
  0,                         /* tp_alloc */
  SparsityInfo_new,                 /* tp_new */

};

// static PyObject* 
// cmeshToolsCommInit(PyObject* self, PyObject* args, PyObject *kwds)
// {
//   int argc;
//   char **argv;
//   PyObject *sys_argv;
//   char *petscDatabaseFilename(0);
//   static char *kwlist[] = {"argv","petscDatabaseFilename"};
//   if(!PyArg_ParseTupleAndKeywords(args,kwds,
// 				  "O|s",kwlist,
// 				  &sys_argv,
// 				  &petscDatabaseFilename))
//     return NULL;
//   argc = PyList_Size(sys_argv);
//   argv = new char* [argc+1];//don't  know why needs one past end cek
//   for (int i=0;i<argc;i++)
//     {
//       argv[i] = PyString_AsString(PyList_GetItem(sys_argv,i));
//     }
//   argv[argc] = new char[1];
//   argv[argc] = '\0';
//   Daetk::Petsc::cc::PetscInitialize(&argc,
// 				    &argv,
// 				    petscDatabaseFilename,
// 				    (char*)("Initializing petsc for Proteus, with options database\n"));
//   PROTEUS_COMM_WORLD = Daetk::Petsc::cc::PETSC_COMM_WORLD;
//   delete [] argv;
//   return Py_None;
// }
// static PyObject* 
// cmeshToolsCommDestroy(PyObject* self, PyObject* args)
// {
//   //mwf add check to make sure not finalized already?
//   PetscTruth finalizedAlready(PETSC_FALSE);
//   int ierr = PetscFinalized(&finalizedAlready);
//   if (!finalizedAlready)
//    PetscFinalize();
//   return Py_None;
// }

}


static PyMethodDef cmeshToolsMethods[] = {
//   {"commInit",
//    (PyCFunction)cmeshToolsCommInit,
//    METH_VARARGS | METH_KEYWORDS,
//    "close petsc"},
//   {"commDestroy",
//     cmeshToolsCommDestroy,
//     METH_VARARGS,
//    "close petsc"},
//   { "globalSum",
//     cmeshToolsGlobalSum,
//     METH_VARARGS, 
//     "sum the value over all subdomains(processes)"},
//   { "globalMax",
//     cmeshToolsGlobalMax,
//     METH_VARARGS, 
//     "take the max of the value over all subdomains(processes)"},
  { "buildPythonMeshInterface",
    cmeshToolsBuildPythonMeshInterface,
    METH_VARARGS, 
    "Provide handles to the C storage of the mesh without any storage (used for global mesh in parallel)"},
  { "buildPythonMeshInterfaceNoArrays",
    cmeshToolsBuildPythonMeshInterfaceNoArrays,
    METH_VARARGS, 
    "Provide handles to the C storage of the mesh without any storage (used for global mesh in parallel)"},
  { "buildPythonMultilevelMeshInterface",
    cmeshToolsBuildPythonMultilevelMeshInterface,
    METH_VARARGS, 
    "Provide handles to the C storage of the multilevelMesh"},
  { "locallyRefineMultilevelMesh",
    cmeshToolsLocallyRefineMultilevelMesh,
    METH_VARARGS, 
    "locally refine mesh given a tag array for elements on current finest level"},
  { "setNewestNodeBases",
    cmeshToolsSetNewestNodeBasesToLongestEdge,
    METH_VARARGS, 
    "set newest node bases to longest edge in case don't have them initialized already"},
  { "generateTetrahedralMeshFromRectangularGrid",
    cmeshToolsGenerateTetrahedralMeshFromRectangularGrid,
    METH_VARARGS, 
    "Build a tetrahedral mesh directly from the description of a rectangular grid"},
  { "computeGeometricInfo_tetrahedron",
    cmeshToolsComputeGeometricInfo_tetrahedron,
    METH_VARARGS, 
    "Compute h, etc."},
  { "allocateGeometricInfo_tetrahedron",
    cmeshToolsAllocateGeometricInfo_tetrahedron,
    METH_VARARGS, 
    "Allocate h, etc."},
  { "generateTriangularMeshFromRectangularGrid",
    cmeshToolsGenerateTriangularMeshFromRectangularGrid,
    METH_VARARGS, 
    "Build a triangular mesh directly from the description of a rectangular grid"},
  { "computeGeometricInfo_triangle",
    cmeshToolsComputeGeometricInfo_triangle,
    METH_VARARGS, 
    "Compute h, etc."},
  { "allocateGeometricInfo_triangle",
    cmeshToolsAllocateGeometricInfo_triangle,
    METH_VARARGS, 
    "Allocate h, etc."},
  { "generateEdgeMeshFromRectangularGrid",
    cmeshToolsGenerateEdgeMeshFromRectangularGrid,
    METH_VARARGS, 
    "Build a edge mesh directly from the description of a rectangular grid"},
  { "computeGeometricInfo_edge",
    cmeshToolsComputeGeometricInfo_edge,
    METH_VARARGS, 
    "Compute h, etc."},
  { "allocateGeometricInfo_edge",
    cmeshToolsAllocateGeometricInfo_edge,
    METH_VARARGS, 
    "Allocate h, etc."},
//   { "partitionElements",
//     cmeshToolsPartitionElements,
//     METH_VARARGS, 
//     "partition the mesh using an element-based partitioning"},
   {"generateFromTriangleMesh",            
   cmeshToolsGenerateFromTriangleMesh,       
   METH_VARARGS,                        
   "just convert from triangle mesh directly"},  /*doc string for method*/
   {"generateFromTriangleFiles",            
   cmeshToolsGenerateFromTriangleFiles,       
   METH_VARARGS,                        
   "just read from triangle node and element files directly"},  /*doc string for method*/
   {"writeTriangleFiles",            
   cmeshToolsWriteTriangleFiles,       
   METH_VARARGS,                        
   "just write out triangle node and element files directly"},  /*doc string for method*/
   {"generateFromTetgenFiles",            
   cmeshToolsGenerateFromTetgenFiles,       
   METH_VARARGS,                        
   "just read from tetgen node and element files directly"},  /*doc string for method*/
   {"generateFromTetgenFilesParallel",            
   cmeshToolsGenerateFromTetgenFilesParallel,       
   METH_VARARGS,                        
   "just read from tetgen node and element files directly, skip some connectivity because it gets rebuilt on subdomains for parallel"},  /*doc string for method*/
   {"writeTetgenFiles",            
   cmeshToolsWriteTetgenFiles,       
   METH_VARARGS,                        
   "just write out tetgen node and element files directly"},  /*doc string for method*/
   {"write3dmFiles",            
   cmeshToolsWrite3dmFiles,       
   METH_VARARGS,                        
   "just write out 3dm node and element files directly"},  /*doc string for method*/
   {"write2dmFiles",            
   cmeshToolsWrite2dmFiles,       
   METH_VARARGS,                        
   "just write out 2dm node and element files directly"},  /*doc string for method*/
   {"generateFrom3DMFile",            
   cmeshToolsGenerateFrom3DMFile,       
   METH_VARARGS,                        
   "just read from 3DM files directly"},  /*doc string for method*/
   {"generateFrom2DMFile",            
   cmeshToolsGenerateFrom2DMFile,       
   METH_VARARGS,                        
   "just read from 2DM files directly"},  /*doc string for method*/
   {"generateFromHexFile",            
   cmeshToolsGenerateFromHexFile,       
   METH_VARARGS,                        
   "just read from Hex files directly"},  /*doc string for method*/
   {"allocateNodeAndElementNodeDataStructures",            
   cmeshToolsAllocateNodeAndElementNodeDataStructures,       
   METH_VARARGS,                        
   "allocate nodeArray, elementNodesArray, elementMaterialTypes, elementNodeMaterialTypes arrays"},  /*doc string for method*/
   {"constructElementBoundaryElementsArray",            
   cmeshToolsConstructElementBoundaryElementsArray,       
   METH_VARARGS,                        
   "generate element boundary information given node and element-> node connectivity"},  /*doc string for method*/
  { "buildLevel0PythonMeshInterface",
    cmeshToolsBuildLevel0PythonMeshInterface,
    METH_VARARGS, 
    "Provide handles to the simplest 'level 0' mesh representation in C "},   
  { "generateHexahedralMeshFromRectangularGrid",
    cmeshToolsGenerateHexahedralMeshFromRectangularGrid,
    METH_VARARGS, 
    "Generates a structured hexahedron"},            
  { "computeGeometricInfo_hexahedron",
    cmeshToolsComputeGeometricInfo_hexahedron,
    METH_VARARGS, 
    "Compute h, etc."},
  { "allocateGeometricInfo_hexahedron",
    cmeshToolsAllocateGeometricInfo_hexahedron,
    METH_VARARGS, 
    "Allocate h, etc."},
  { "generateNURBSMeshFromRectangularGrid",
    cmeshToolsGenerateNURBSMeshFromRectangularGrid,
    METH_VARARGS, 
    "Generates a structured NURBS"},  
  { "computeGeometricInfo_NURBS",
    cmeshToolsComputeGeometricInfo_NURBS,
    METH_VARARGS, 
    "Compute h, etc."},
  { "allocateGeometricInfo_NURBS",
    cmeshToolsAllocateGeometricInfo_NURBS,
    METH_VARARGS, 
    "Allocate h, etc."},
  { "deleteMeshDataStructures",
    cmeshToolsDeleteMeshDataStructures,
    METH_VARARGS, 
    "Allocate h, etc."},
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcmeshTools(void)
{
  PyObject *m,*d;
  if (PyType_Ready(&CMeshType) < 0)
    return;
//   if (PyType_Ready(&ParVecType) < 0)
//     return;
//   if (PyType_Ready(&ParMatType) < 0)
//     return;
//   if (PyType_Ready(&CKSPType) < 0)
//     return;
  if (PyType_Ready(&CMultilevelMeshType) < 0)
    return;
  if (PyType_Ready(&SparsityInfoType) < 0)
    return;
  m = Py_InitModule3("cmeshTools", 
                     cmeshToolsMethods,
                     "cmeshTools module");
//   Py_INCREF(&CKSPType);
//   PyModule_AddObject(m, "KSP", (PyObject *)&CKSPType);
//   Py_INCREF(&ParVecType);
//   PyModule_AddObject(m, "ParVec", (PyObject *)&ParVecType);
//   Py_INCREF(&ParMatType);
//   PyModule_AddObject(m, "ParMat", (PyObject *)&ParMatType);
  Py_INCREF(&CMeshType);
  PyModule_AddObject(m, "CMesh", (PyObject *)&CMeshType);
  Py_INCREF(&CMultilevelMeshType);
  PyModule_AddObject(m, "CMultilevelMesh", (PyObject *)&CMultilevelMeshType);
  Py_INCREF(&SparsityInfoType);
  PyModule_AddObject(m, "SparsityInfo", (PyObject *)&SparsityInfoType);
  d = PyModule_GetDict(m);
  import_array();
//   import_flcbdfWrappers();
}
}
/** @} */
