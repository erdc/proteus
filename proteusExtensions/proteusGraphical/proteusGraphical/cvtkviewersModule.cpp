#include "Python.h"
#include "numpy/arrayobject.h"
#include "vtkPythonUtil.h"
#include "vtkviewers.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "mesh.h"
#include "cmeshToolsModule.h"

#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

#ifdef USE_VTK_5_6_OR_EARLIER
# define VTK_PYTHON_OBJECT_RETURN(X) vtkPythonGetObjectFromPointer(X)
#else
# define VTK_PYTHON_OBJECT_RETURN(X) vtkPythonUtil::GetObjectFromPointer(X)
#endif

#ifdef USE_VTK_5_9_OR_LATER
#define VTK_PYTHON_POINTER_RETURN(X,Y) vtkPythonUtil::GetPointerFromObject(X,Y)
#else
#define VTK_PYTHON_POINTER_RETURN(X,Y) vtkPythonGetPointerFromObject(X,Y)
#endif

extern "C"
{
  static PyObject* 
  cvtkviewersGetUnstructuredGridFromMesh(PyObject* self,
					 PyObject* args)
  {
    PyObject * nodeArray, *elementNodesArray, *elementMaterialTypes = 0, *nodeMaterialTypes = 0; 
    vtkUnstructuredGrid* vtkUgrid(0);
    
    if (!PyArg_ParseTuple(args,
			  "OO|OO",
			  &nodeArray,
			  &elementNodesArray,
			  &elementMaterialTypes,
			  &nodeMaterialTypes))
      return NULL;
    if (nodeMaterialTypes != 0)
      assert(elementMaterialTypes != 0);
    if (elementMaterialTypes != 0)
      {
	if (nodeMaterialTypes != 0)
	  {
	    vtkUgrid = vtkUnstructuredGridFromMesh(SHAPE(nodeArray)[0],
						   SHAPE(elementNodesArray)[0],
						   SHAPE(elementNodesArray)[1],
						   DDATA(nodeArray),
						   IDATA(elementNodesArray),
						   IDATA(elementMaterialTypes),
						   IDATA(nodeMaterialTypes));
	  }
	else
	  {
	    vtkUgrid = vtkUnstructuredGridFromMesh(SHAPE(nodeArray)[0],
						   SHAPE(elementNodesArray)[0],
						   SHAPE(elementNodesArray)[1],
						   DDATA(nodeArray),
						   IDATA(elementNodesArray),
						   IDATA(elementMaterialTypes));
	  }
      }
    else
      {
	vtkUgrid = vtkUnstructuredGridFromMesh(SHAPE(nodeArray)[0],
					       SHAPE(elementNodesArray)[0],
					       SHAPE(elementNodesArray)[1],
					       DDATA(nodeArray),
					       IDATA(elementNodesArray));
      }
    return VTK_PYTHON_OBJECT_RETURN(vtkUgrid);
  }

  static PyObject* 
  cvtkviewersGetUnstructuredGridFromQuadraticTriangleMesh(PyObject* self,
							  PyObject* args)
  {
    int numNodes, numNodesTotal, nElem;
    PyObject *nodeArray, *l2g, *edgeNodesArray; 
    
    if (!PyArg_ParseTuple(args,
			  "iiiOOO",
			  &numNodes,
			  &numNodesTotal,
			  &nElem,
			  &nodeArray,
			  &l2g,
			  &edgeNodesArray))
      {
	return NULL;
      }
    
    vtkUnstructuredGrid* vtkUgrid = vtkUnstructuredGridFromQuadraticTriangleMesh(numNodes,
										 numNodesTotal,
										 nElem,
										 DDATA(nodeArray),
										 IDATA(l2g),
										 IDATA(edgeNodesArray));
    
    return VTK_PYTHON_OBJECT_RETURN(vtkUgrid);
  }
  
  static PyObject* 
  cvtkviewersGetUnstructuredGridFromQuadraticTetMesh(PyObject* self,
						     PyObject* args)
  {
    int numNodes, numNodesTotal, nElem;
    PyObject *nodeArray, *l2g, *edgeNodesArray; 
    
    if (!PyArg_ParseTuple(args,
			  "iiiOOO",
			  &numNodes,
			  &numNodesTotal,
			  &nElem,
			  &nodeArray,
			  &l2g,
			  &edgeNodesArray))
      {
	return NULL;
      }
    
    vtkUnstructuredGrid* vtkUgrid = vtkUnstructuredGridFromQuadraticTetMesh(numNodes,
									    numNodesTotal,
									    nElem,
									    DDATA(nodeArray),
									    IDATA(l2g),
									    IDATA(edgeNodesArray));
    
    return VTK_PYTHON_OBJECT_RETURN(vtkUgrid);
  }
  
  static PyObject* 
  cvtkviewersPrepareScalarValueArray(PyObject* self, PyObject* args)
  {
    int length=1;
    PyObject *scalarsArray;
    
    if (!PyArg_ParseTuple(args,
			  "O",
			  &scalarsArray))
      return NULL;
    for(int i=0;i<ND(scalarsArray);i++)
      length*=SHAPE(scalarsArray)[i];
    vtkDoubleArray* vtkScalars = vtkPrepareScalarValueArray(length,
							    DDATA(scalarsArray));
    
    return VTK_PYTHON_OBJECT_RETURN(vtkScalars);
  }
  
  static PyObject* 
  cvtkviewersPrepareVectorValueArray(PyObject* self, PyObject* args)
  {
    PyObject *vectorsArray;
    
    if (!PyArg_ParseTuple(args,
			  "O",
			  &vectorsArray))
      return NULL;
    vtkDoubleArray* vtkVectors = vtkPrepareVectorValueArray(SHAPE(vectorsArray)[0],
    							    3,
    							    DDATA(vectorsArray));
    
    return VTK_PYTHON_OBJECT_RETURN(vtkVectors);
  }
  
  
  static PyObject* 
  cvtkviewersPrepareVTKPoints(PyObject* self, PyObject* args)
  {
    PyObject *nodeArray;
    
    if (!PyArg_ParseTuple(args,
			  "O",
			  &nodeArray))
      return NULL;
    
    vtkPoints* points= vtkPrepareVTKPoints(SHAPE(nodeArray)[0],
					   DDATA(nodeArray));
    
    return VTK_PYTHON_OBJECT_RETURN(points);
  }
  
  static PyObject* 
  cvtkviewersPrepareVTKPointsFrom1DArray(PyObject* self, PyObject* args)
  {
    PyObject *nodeArray;
    
    if (!PyArg_ParseTuple(args,
			  "O",
			  &nodeArray))
      return NULL;
    
    vtkPoints* points= vtkPrepareVTKPoints(SHAPE(nodeArray)[0]/3,
					   DDATA(nodeArray));
    
    return VTK_PYTHON_OBJECT_RETURN(points);
  }
  
  static PyObject* 
  cvtkviewersPrepareVTKPoints3(PyObject* self, PyObject* args)
  {
    int length=1;
    PyObject *nodeArray;
    
    if (!PyArg_ParseTuple(args,
			  "O",
			  &nodeArray))
      return NULL;
    for(int i=0;i<(ND(nodeArray)-1);i++)
      length*=SHAPE(nodeArray)[i];
    vtkPoints* points= vtkPrepareVTKPoints(length,
					   DDATA(nodeArray));
    
    return VTK_PYTHON_OBJECT_RETURN(points);
  }
    
  static PyObject* 
  cvtkviewersGetPolyDataBoundaryMesh(PyObject* self,
				     PyObject* args)
  {
    PyObject * nodeArray, *elementBoundaryNodesArray,
      * elementBoundaryMaterialTypes; 
    
    if (!PyArg_ParseTuple(args,
			  "OOO",
			  &nodeArray,
			  &elementBoundaryNodesArray,
			  &elementBoundaryMaterialTypes))
      return NULL;
    vtkPolyData* vtkBgrid = vtkPolyDataBoundaryMesh(SHAPE(nodeArray)[0],
						    SHAPE(elementBoundaryNodesArray)[0],
						    SHAPE(elementBoundaryNodesArray)[1],
						    DDATA(nodeArray),
						    IDATA(elementBoundaryNodesArray),
						    IDATA(elementBoundaryMaterialTypes));
    
    return VTK_PYTHON_OBJECT_RETURN(vtkBgrid);
  }
  
  static PyObject* 
  cvtkviewersGetMeshFromVTKUnstructuredGrid(PyObject* self,
					    PyObject* args)
  {
    PyObject *vtkMeshIn, * cmesh;
    int copyIdsFlag = 0;
    if (!PyArg_ParseTuple(args,
			  "OO|i",
			  &vtkMeshIn,
			  &cmesh,
			  &copyIdsFlag))
      return NULL;
    vtkUnstructuredGrid * vtkMesh = (vtkUnstructuredGrid*) VTK_PYTHON_POINTER_RETURN(vtkMeshIn,"vtkUnstructuredGrid");
    assert(vtkMesh);
    bool copyMaterialIds = (copyIdsFlag == 1); 
    int defaultElementMaterialFlag(0),defaultNodeMaterialFlag(0);
    
    bool failed = meshElementAndNodeArraysFromVTKUnstructuredGrid(vtkMesh,
								  MESH(cmesh).nNodes_global,
								  MESH(cmesh).nElements_global,
								  MESH(cmesh).nNodes_element,
								  MESH(cmesh).nodeArray,
								  MESH(cmesh).elementNodesArray,
								  MESH(cmesh).elementMaterialTypes,
								  MESH(cmesh).nodeMaterialTypes,
								  copyMaterialIds,
								  defaultElementMaterialFlag,
								  defaultNodeMaterialFlag);
    if (failed)
      assert(!failed);
    
    if (MESH(cmesh).nNodes_element == 2)
      constructElementBoundaryElementsArray_edge(MESH(cmesh));
    else if (MESH(cmesh).nNodes_element == 3)
      constructElementBoundaryElementsArray_triangle(MESH(cmesh));
    else if (MESH(cmesh).nNodes_element == 4)
      constructElementBoundaryElementsArray_tetrahedron(MESH(cmesh));
    Py_INCREF(Py_None); 
    return Py_None;
  }
  static PyObject* 
  cvtkviewersClassifyElementMaterialPropertiesFromVTKUnstructuredGridSolid(PyObject* self,
									   PyObject* args)
  {
    PyObject *vtkSolidIn,*elementBarycentersArray,*elementDiametersArray,*elementMaterialTypes;
    int nElements,newMaterialId,verbose = 0;
    double tol_hFactor;
    if (!PyArg_ParseTuple(args,
			  "iOOOOd|i",
			  &newMaterialId,
			  &vtkSolidIn,
			  &elementBarycentersArray,
			  &elementDiametersArray,
			  &elementMaterialTypes,
			  &tol_hFactor,
			  &verbose))
      return NULL;
    vtkUnstructuredGrid * vtkSolid = (vtkUnstructuredGrid*) VTK_PYTHON_POINTER_RETURN(vtkSolidIn,"vtkUnstructuredGrid");
    assert(vtkMesh);
    nElements = SHAPE(elementMaterialTypes)[0];
    bool failed = classifyElementMaterialPropertiesFromVTKUnstructuredGridSolid(vtkSolid,
										nElements,
										newMaterialId,
										DDATA(elementBarycentersArray),
										DDATA(elementDiametersArray),
										IDATA(elementMaterialTypes),
										tol_hFactor,
										verbose);
    if (failed)
      assert(!failed);

    Py_INCREF(Py_None); 
    return Py_None;
  }
  static PyObject* 
  cvtkviewersClassifyElementMaterialPropertiesFromVTKUnstructuredGridNeighborhood(PyObject* self,
										  PyObject* args)
  {
    PyObject *vtkSolidIn,*elementBarycentersArray,*elementDiametersArray,*elementMaterialTypes;
    int nElements,newMaterialId,verbose = 0,useAbsoluteTolerance;
    double tolerance;
    if (!PyArg_ParseTuple(args,
			  "iOOOOid|i",
			  &newMaterialId,
			  &vtkSolidIn,
			  &elementBarycentersArray,
			  &elementDiametersArray,
			  &elementMaterialTypes,
			  &useAbsoluteTolerance,
			  &tolerance,
			  &verbose))
      return NULL;
    vtkUnstructuredGrid * vtkSolid = (vtkUnstructuredGrid*) VTK_PYTHON_POINTER_RETURN(vtkSolidIn,"vtkUnstructuredGrid");
    assert(vtkMesh);
    nElements = SHAPE(elementMaterialTypes)[0];
    bool failed = classifyElementMaterialPropertiesFromVTKUnstructuredGridNeighborhood(vtkSolid,
										       nElements,
										       newMaterialId,
										       DDATA(elementBarycentersArray),
										       DDATA(elementDiametersArray),
										       IDATA(elementMaterialTypes),
										       useAbsoluteTolerance,
										       tolerance,
										       verbose);
    if (failed)
      assert(!failed);

    Py_INCREF(Py_None); 
    return Py_None;
  }
  
  static PyMethodDef cvtkviewersMethods[] = {
    
    {"getUnstructuredGridFromMesh",
     cvtkviewersGetUnstructuredGridFromMesh,
     METH_VARARGS,
     "generate vtkUnstructuredGrid data structure directly in c++"},
    {"getUnstructuredGridFromQuadraticTriangleMesh",
     cvtkviewersGetUnstructuredGridFromQuadraticTriangleMesh,
     METH_VARARGS,
     "generate vtkUnstructuredGrid data structure directly in c++"},
    {"getUnstructuredGridFromQuadraticTetMesh",
     cvtkviewersGetUnstructuredGridFromQuadraticTetMesh,
     METH_VARARGS,
     "generate vtkUnstructuredGrid data structure directly in c++"},
    {"prepareScalarValueArray",
     cvtkviewersPrepareScalarValueArray,
     METH_VARARGS,
     "generate vtkDoubleArray data structure full of 1D touples directly in c++"},
    {"prepareVectorValueArray",
     cvtkviewersPrepareVectorValueArray,
     METH_VARARGS,
     "generate vtkDoubleArray data structure full of multi-dimension touples directly in c++"},
    {"prepareVTKPointsFrom1DArray",
     cvtkviewersPrepareVTKPointsFrom1DArray,
     METH_VARARGS,
     "generate ctkPoints data structure directly in c++"},
    {"prepareVTKPoints",
     cvtkviewersPrepareVTKPoints,
     METH_VARARGS,
     "generate ctkPoints data structure directly in c++"},
    {"prepareVTKPoints3",
     cvtkviewersPrepareVTKPoints3,
     METH_VARARGS,
     "generate ctkPoints data structure directly in c++ where the number of points = len(array)/3"},
    {"getPolyDataBoundaryMesh",
     cvtkviewersGetPolyDataBoundaryMesh,
     METH_VARARGS,
     "generate vtkPolyData representation of boundary mesh in c++"},
    {"getMeshFromVTKUnstructuredGrid",
     cvtkviewersGetMeshFromVTKUnstructuredGrid,
     METH_VARARGS,
     "generate meshTools mesh representation from vtkUnstructured grid"},
    {"classifyElementMaterialPropertiesFromVTKUnstructuredGridSolid",
     cvtkviewersClassifyElementMaterialPropertiesFromVTKUnstructuredGridSolid,
     METH_VARARGS,
     "assign new material id to elements in mesh if barycenter in (or close?) to solid defined by vtkUnstructuredGrid"},
     {"classifyElementMaterialPropertiesFromVTKUnstructuredGridNeighborhood",
     cvtkviewersClassifyElementMaterialPropertiesFromVTKUnstructuredGridNeighborhood,
     METH_VARARGS,
     "assign new material id to elements in mesh if barycenter within spherical neighborhood of solid defined by vtkUnstructuredGrid"},
    { NULL,NULL,0,NULL}
  };
  
  PyMODINIT_FUNC initcvtkviewers(void)
  {
    PyObject *m,*d;
    m = Py_InitModule("cvtkviewers", cvtkviewersMethods);
    d = PyModule_GetDict(m);
    import_array();
  }
}//extern "C"
