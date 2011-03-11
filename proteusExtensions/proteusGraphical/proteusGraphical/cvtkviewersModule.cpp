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

extern "C"
{
  static PyObject* 
  cvtkviewersGetUnstructuredGridFromMesh(PyObject* self,
					 PyObject* args)
  {
    PyObject * nodeArray, *elementNodesArray, *elementMaterialTypes = 0; 
    vtkUnstructuredGrid* vtkUgrid(0);
    
    if (!PyArg_ParseTuple(args,
			  "OO|O",
			  &nodeArray,
			  &elementNodesArray,
			  &elementMaterialTypes))
      return NULL;
    if (elementMaterialTypes != 0)
      {
	vtkUgrid = vtkUnstructuredGridFromMesh(SHAPE(nodeArray)[0],
					       SHAPE(elementNodesArray)[0],
					       SHAPE(elementNodesArray)[1],
					       DDATA(nodeArray),
					       IDATA(elementNodesArray),
					       IDATA(elementMaterialTypes));
      }
    else
      {
	vtkUgrid = vtkUnstructuredGridFromMesh(SHAPE(nodeArray)[0],
					       SHAPE(elementNodesArray)[0],
					       SHAPE(elementNodesArray)[1],
					       DDATA(nodeArray),
					       IDATA(elementNodesArray));
      }
    return vtkPythonUtil::GetObjectFromPointer(vtkUgrid);
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
    
    return vtkPythonUtil::GetObjectFromPointer(vtkUgrid);
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
    
    return vtkPythonUtil::GetObjectFromPointer(vtkUgrid);
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
    
    return vtkPythonUtil::GetObjectFromPointer(vtkScalars);
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
    
    return vtkPythonUtil::GetObjectFromPointer(vtkVectors);
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
    
    return vtkPythonUtil::GetObjectFromPointer(points);
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
    
    return vtkPythonUtil::GetObjectFromPointer(points);
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
    
    return vtkPythonUtil::GetObjectFromPointer(points);
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
    
    return vtkPythonUtil::GetObjectFromPointer(vtkBgrid);
  }
  
  static PyObject* 
  cvtkviewersGetMeshFromVTKUnstructuredGrid(PyObject* self,
					    PyObject* args)
  {
    PyObject *vtkMeshIn, * cmesh;
    
    if (!PyArg_ParseTuple(args,
			  "OO",
			  &vtkMeshIn,
			  &cmesh))
      return NULL;
    vtkUnstructuredGrid * vtkMesh = (vtkUnstructuredGrid*) vtkPythonUtil::GetPointerFromObject(vtkMeshIn,"vtkUnstructuredGrid");
    assert(vtkMesh);
    bool copyMaterialIds = false; 
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
