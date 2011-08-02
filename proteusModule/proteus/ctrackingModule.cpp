#include "Python.h"
#include "numpy/arrayobject.h"
#include "tracking.h"
#include "superluWrappersModule.h"
#include <iostream>

#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define CSRVAL(p) ((double*)((SparseMatrix*)p)->A.nzval)

#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define CSRVAL(p) ((double*)((SparseMatrix*)p)->A.nzval)

static PyObject* ctracking_trackPointsConstVelocity1d(PyObject* self,
						      PyObject* args)
{
  int nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nPointsToTrack;

  double cVelocity,
    dir,
    zeroTolForTracking;

  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *x_depart_times,
    *x_arrive_times,
    *x_in,
    *x_element,
    *x_out,
    *flag;

  if (!PyArg_ParseTuple(args,
			"iiiiOOOddOOidOOOO",
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&cVelocity,             
			&dir,                   
			&x_depart_times,               
			&x_arrive_times,               
			&nPointsToTrack,        
			&zeroTolForTracking,    
			&x_in,                
			&x_element,           
			&x_out,               
			&flag))    
    return NULL;

  trackPointsConstVelocity1d(nElements_global,               //mesh representation
			     nNodes_global,
			     nNodes_element,
			     nElementBoundaries_element,
			     DDATA(nodeArray),
			     IDATA(elementNodesArray),
			     IDATA(elementNeighborsArray), //local boundary id is associated with node across from boundary 
			     cVelocity,                  //characteristic speed (velocity) representation
			     dir,                        //direction in time
			     DDATA(x_depart_times),      //in  -- time of departure for each point
                             DDATA(x_arrive_times),      //desired stopping time
   		                                         //out -- stopping time for each point 
			     nPointsToTrack,             //number of points to track
			     zeroTolForTracking,         //ignore point if |u| < eps or |v| < eps
			     DDATA(x_in),                //points for tracking (always 3d)
			     IDATA(x_element),           //in -- element where point i is located at tIn
				                         //out -- element where point i is located at tOut
			     DDATA(x_out),               //stopping location for point i at tOut
			     IDATA(flag));               //in: > -2  -- track point
                                                         //in: -3    -- skip point
                                                         //out: -1   -- point in interior at tOut
					                 //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                         //out: -3   -- did not track (e.g.,  or u = 0)

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctracking_trackPointsC0P1Velocity1d(PyObject* self,
						   PyObject* args)
{
  int nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nPointsToTrack;

  double dir,
    zeroTolForTracking;

  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *cvelocity_l2g,
    *cvelocity_dof,
    *x_depart_times,
    *x_arrive_times,
    *x_in,
    *x_element,
    *x_out,
    *flag;

  if (!PyArg_ParseTuple(args,
			"iiiiOOOOOdOOidOOOO",
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&cvelocity_l2g,
			&cvelocity_dof,
			&dir,                   
			&x_depart_times,               
			&x_arrive_times,               
			&nPointsToTrack,        
			&zeroTolForTracking,    
			&x_in,                
			&x_element,           
			&x_out,               
			&flag))    
    return NULL;

  trackPointsC0P1Velocity1d(nElements_global,               //mesh representation
			    nNodes_global,
			    nNodes_element,
			    nElementBoundaries_element,
			    DDATA(nodeArray),
			    IDATA(elementNodesArray),
			    IDATA(elementNeighborsArray), //local boundary id is associated with node across from boundary 
			    IDATA(cvelocity_l2g),
			    DDATA(cvelocity_dof),       //characteristic speed (velocity) representation
			    dir,                        //direction in time
			    DDATA(x_depart_times),      //in  -- time of departure for each point
			    DDATA(x_arrive_times),      //desired stopping time
			    //out -- stopping time for each point 
			    nPointsToTrack,             //number of points to track
			    zeroTolForTracking,         //ignore point if |u| < eps or |v| < eps
			    DDATA(x_in),                //points for tracking (always 3d)
			    IDATA(x_element),           //in -- element where point i is located at tIn
			    //out -- element where point i is located at tOut
			    DDATA(x_out),               //stopping location for point i at tOut
			    IDATA(flag));                        //in: > -2  -- track point
                                                                 //in: -3    -- skip point
                                                                 //out: -1   -- point in interior at tOut
  			                                         //out: -2   -- point exited domain somewhere in (tIn,tOut)
                                                                 //out: -3   -- did not track (e.g.,  or u = 0)
 
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctracking_trackPointsRT0Velocity2d(PyObject* self,
						    PyObject* args)
{
  int localVelocityRepresentationFlag,
    nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nPointsToTrack,
    debugLevel;
  
  double dir,
    zeroTolForTracking;

  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *elementBoundariesArray,
    *elementBoundaryBarycentersArray,
    *elementLocalBoundaryOuterNormalsArray,
    *cvelocity_l2g,
    *cvelocity_dof,
    *x_depart_times,
    *x_arrive_times,
    *x_in,
    *x_element,
    *x_out,
    *flag;

  debugLevel = 0;
  if (!PyArg_ParseTuple(args,
			"iiiiiOOOOOOOOdOOidOOOO|i",
			&localVelocityRepresentationFlag,
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&elementBoundariesArray,
			&elementBoundaryBarycentersArray,
			&elementLocalBoundaryOuterNormalsArray,
			&cvelocity_l2g,
			&cvelocity_dof,
			&dir,                   
			&x_depart_times,               
			&x_arrive_times,               
			&nPointsToTrack,        
			&zeroTolForTracking,    
			&x_in,                
			&x_element,           
			&x_out,               
			&flag,
			&debugLevel))    
    return NULL;

  trackPointsRT0Velocity2d(debugLevel,
			   localVelocityRepresentationFlag,
			   nElements_global,             
			   nNodes_global,
			   nNodes_element,
			   nElementBoundaries_element,
			   DDATA(nodeArray),
			   IDATA(elementNodesArray),
			   IDATA(elementNeighborsArray), 
			   IDATA(elementBoundariesArray),
			   DDATA(elementBoundaryBarycentersArray),
			   DDATA(elementLocalBoundaryOuterNormalsArray),
			   IDATA(cvelocity_l2g),
			   DDATA(cvelocity_dof),       
			   dir,                        
			   DDATA(x_depart_times),      
			   DDATA(x_arrive_times),      
			   nPointsToTrack,             
			   zeroTolForTracking,         
			   DDATA(x_in),                
			   IDATA(x_element),           
			   DDATA(x_out),
			   IDATA(flag));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* ctracking_trackPointsRT0Velocity2dWithTrajectories(PyObject* self,
								    PyObject* args)
{
  int localVelocityRepresentationFlag,
    nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nPointsToTrack,
    debugLevel;
  
  double dir,
    zeroTolForTracking;

  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *elementBoundariesArray,
    *elementBoundaryBarycentersArray,
    *elementLocalBoundaryOuterNormalsArray,
    *cvelocity_l2g,
    *cvelocity_dof,
    *x_depart_times,
    *x_arrive_times,
    *x_in,
    *x_element,
    *x_out,
    *flag;

  //try to build arrays for trajectories and return
  int n_traj,n_tracked;
  double * x_traj = 0, * t_traj = 0;
  int * e_traj = 0, * offsets_traj = 0;
  debugLevel = 0;
  if (!PyArg_ParseTuple(args,
			"iiiiiOOOOOOOOdOOidOOOO|i",
			&localVelocityRepresentationFlag,
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&elementBoundariesArray,
			&elementBoundaryBarycentersArray,
			&elementLocalBoundaryOuterNormalsArray,
			&cvelocity_l2g,
			&cvelocity_dof,
			&dir,                   
			&x_depart_times,               
			&x_arrive_times,               
			&nPointsToTrack,        
			&zeroTolForTracking,    
			&x_in,                
			&x_element,           
			&x_out,               
			&flag,
			&debugLevel))    
    return NULL;

  trackPointsRT0Velocity2dWithTrajectories(debugLevel,
					   localVelocityRepresentationFlag,
					   nElements_global,             
					   nNodes_global,
					   nNodes_element,
					   nElementBoundaries_element,
					   DDATA(nodeArray),
					   IDATA(elementNodesArray),
					   IDATA(elementNeighborsArray), 
					   IDATA(elementBoundariesArray),
					   DDATA(elementBoundaryBarycentersArray),
					   DDATA(elementLocalBoundaryOuterNormalsArray),
					   IDATA(cvelocity_l2g),
					   DDATA(cvelocity_dof),       
					   dir,                        
					   DDATA(x_depart_times),      
					   DDATA(x_arrive_times),      
					   nPointsToTrack,             
					   zeroTolForTracking,         
					   DDATA(x_in),                
					   IDATA(x_element),           
					   DDATA(x_out),
					   IDATA(flag),
					   n_traj,
					   n_tracked,
					   offsets_traj,
					   x_traj,
					   t_traj,
					   e_traj);
  npy_intp dim[1],dim2[2]; 
  dim[0] = n_traj;
  PyArrayObject *t_traj_py = (PyArrayObject *)PyArray_SimpleNew(1,dim,PyArray_DOUBLE);
  double* t_ptr = DDATA(t_traj_py);

  dim2[0] = n_traj; dim2[1] = 3;
  PyArrayObject *x_traj_py = (PyArrayObject *)PyArray_SimpleNew(2,dim2,PyArray_DOUBLE);
  double* x_ptr = DDATA(x_traj_py);

  dim[0] = n_traj;
  PyArrayObject *e_traj_py = (PyArrayObject *)PyArray_SimpleNew(1,dim,PyArray_INT);
  int* e_ptr = IDATA(e_traj_py);

  dim[0] = n_tracked+1;
  PyArrayObject *o_traj_py = (PyArrayObject *)PyArray_SimpleNew(1,dim,PyArray_INT);
  int* o_ptr = IDATA(o_traj_py);

  for (int i=0; i < n_traj; i++)
    {
      t_ptr[i] = t_traj[i]; e_ptr[i] = e_traj[i]; 
      for (int j=0; j < 3; j++)
	x_ptr[i*3+j] = x_traj[i*3+j];
    }
  for (int i=0; i < n_tracked+1; i++)
    {
      o_ptr[i] = offsets_traj[i];
    }
  return Py_BuildValue("(O,O,O,O)",PyArray_Return(x_traj_py),PyArray_Return(t_traj_py),PyArray_Return(e_traj_py),PyArray_Return(o_traj_py));

}

static PyObject* ctracking_getOuterNormals_affineSimplex(PyObject* self,
							 PyObject* args)
{
  int nSpace;
  


  PyObject *boundaryNormals,
    *jacobianInverseArray,
    *unitNormalArray;

  if (!PyArg_ParseTuple(args,
			"OOO",
			&boundaryNormals,
			&jacobianInverseArray,
			&unitNormalArray))
    return NULL;
  nSpace = SHAPE(boundaryNormals)[0] -1;
  if (nSpace == 1)
    {
      getOuterNormals_affineSimplex_1d(SHAPE(jacobianInverseArray)[0],
				       SHAPE(unitNormalArray)[1],
				       SHAPE(jacobianInverseArray)[1],
				       DDATA(boundaryNormals),
				       DDATA(jacobianInverseArray),
				       DDATA(unitNormalArray));

    }
  else if (nSpace == 2)
    {
      getOuterNormals_affineSimplex_2d(SHAPE(jacobianInverseArray)[0],
				       SHAPE(unitNormalArray)[1],
				       SHAPE(jacobianInverseArray)[1],
				       DDATA(boundaryNormals),
				       DDATA(jacobianInverseArray),
				       DDATA(unitNormalArray));

    }
  else
    {
      assert(nSpace == 3);
      getOuterNormals_affineSimplex_3d(SHAPE(jacobianInverseArray)[0],
				       SHAPE(unitNormalArray)[1],
				       SHAPE(jacobianInverseArray)[1],
				       DDATA(boundaryNormals),
				       DDATA(jacobianInverseArray),
				       DDATA(unitNormalArray));

    }
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* ctracking_setNodeOnBoundaryArray(PyObject* self,
						  PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryNodesArray,
    *nodeOnBoundaryArray;
  if (!PyArg_ParseTuple(args,
			"OOO",
			&exteriorElementBoundariesArray,
			&elementBoundaryNodesArray,
			&nodeOnBoundaryArray))
    return NULL;

  setNodeOnBoundaryArray(SHAPE(exteriorElementBoundariesArray)[0],
			 SHAPE(elementBoundaryNodesArray)[1],
			 IDATA(exteriorElementBoundariesArray),
			 IDATA(elementBoundaryNodesArray),
			 IDATA(nodeOnBoundaryArray));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef ctrackingMethods[] = {
 { "trackPointsConstVelocity1d",
   ctracking_trackPointsConstVelocity1d,
   METH_VARARGS, 
   "track points on unstructured 1d mesh assuming constant velocity"},
 { "trackPointsC0P1Velocity1d",
   ctracking_trackPointsC0P1Velocity1d,
   METH_VARARGS, 
   "track points on unstructured 1d mesh assuming C0 P1 velocity"},
 { "trackPointsRT0Velocity2d",
   ctracking_trackPointsRT0Velocity2d,
   METH_VARARGS, 
   "track points on unstructured 1d mesh assuming RT0 velocity"},
 { "trackPointsRT0Velocity2dWithTrajectories",
   ctracking_trackPointsRT0Velocity2dWithTrajectories,
   METH_VARARGS, 
   "track points on unstructured 2d mesh assuming RT0 velocity, returns particle trajectory history too"},
 { "getOuterNormals_affineSimplex",
   ctracking_getOuterNormals_affineSimplex,
   METH_VARARGS, 
   "grab outer normals for triangular mesh"},
 { "setNodeOnBoundaryArray",
   ctracking_setNodeOnBoundaryArray,
   METH_VARARGS, 
   "tag nodes on exterior boundary"},
 { NULL,NULL,0,NULL}
};

extern "C"
{
PyMODINIT_FUNC initctracking(void)
{
  PyObject *m,*d;
  m = Py_InitModule("ctracking", ctrackingMethods);
  d = PyModule_GetDict(m);
  import_array();
}
}//extern "C"
