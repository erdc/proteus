#include "Python.h"
#include "numpy/arrayobject.h"
#include "timeIntegration.h"
/** \file ctimeIntegrationModule.c
    \defgroup ctimeIntegration ctimeIntegration
    \brief Python interfact to timeIntegration library 
    @{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

static PyObject*
ctimeIntegrationPsiTCtteDT(PyObject* self,PyObject* args)
{
  double tau,dtn,dtnm1,dtnp1;
  PyObject *yn,*ypn,*ypnm1,*rval;
  int nPoints = 1;
  if(!PyArg_ParseTuple(args,"dddOOO",
                       &tau,
		       &dtn,
		       &dtnm1,
		       &yn,
		       &ypn,
		       &ypnm1))
    return NULL;
  nPoints = SHAPE(yn)[0]*SHAPE(yn)[1];

  psiTCtteDT(nPoints,
	     tau,
	     dtn,
	     dtnm1,
	     DDATA(yn),
	     DDATA(ypn),
	     DDATA(ypnm1),
	     &dtnp1);

  /*Py_INCREF(Py_None);*/
  rval = Py_BuildValue("d",dtnp1);
  return rval;

}

static PyObject*
ctimeIntegrationApplyDGlimitingP1Lagrange1d(PyObject* self,
					    PyObject* args)
{
  PyObject *elementNodesArray,*elementNeighborsArray,*nodeArray,
    *elementBarycentersArray,*l2g, *tag, *Uin, *Uout;
  int limiterFlag = 1; /*type of limiting, muscl by default*/
  if(!PyArg_ParseTuple(args,"OOOOOOOO|i",
                       &elementNodesArray,
		       &elementNeighborsArray,
		       &nodeArray,
		       &elementBarycentersArray,
		       &l2g,
		       &tag,
		       &Uin,
		       &Uout,
		       &limiterFlag))
    return NULL;
  
  applyDGlimitingP1Lagrange1d(limiterFlag,
			      SHAPE(elementNodesArray)[0],
			      SHAPE(elementNodesArray)[1],
			      SHAPE(elementNeighborsArray)[1],
			      SHAPE(l2g)[1], 
			      IDATA(elementNodesArray),
			      IDATA(elementNeighborsArray),
			      DDATA(nodeArray),
			      DDATA(elementBarycentersArray),
			      IDATA(l2g),
			      IDATA(tag),
			      DDATA(Uin),
			      DDATA(Uout));

  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject*
ctimeIntegrationApplyDGlimitingP1Lagrange1d_withVacuumTol(PyObject* self,
							  PyObject* args)
{
  PyObject *elementNodesArray,*elementNeighborsArray,*nodeArray,
    *elementBarycentersArray,*l2g, *tag, *Uin, *Uout;
  int enforcePositivity = 0; /*switch to conservative limiting for u < tol or |u| < tol*/
  double vacuumTol  =0.0;
  if(!PyArg_ParseTuple(args,"idOOOOOOOO",
		       &enforcePositivity,
		       &vacuumTol,
                       &elementNodesArray,
		       &elementNeighborsArray,
		       &nodeArray,
		       &elementBarycentersArray,
		       &l2g,
		       &tag,
		       &Uin,
		       &Uout))
    return NULL;
  
  applyDGlimitingP1Lagrange1d_withVacuumTol(enforcePositivity,
					    vacuumTol,
					    SHAPE(elementNodesArray)[0],
					    SHAPE(elementNodesArray)[1],
					    SHAPE(elementNeighborsArray)[1],
					    SHAPE(l2g)[1], 
					    IDATA(elementNodesArray),
					    IDATA(elementNeighborsArray),
					    DDATA(nodeArray),
					    DDATA(elementBarycentersArray),
					    IDATA(l2g),
					    IDATA(tag),
					    DDATA(Uin),
					    DDATA(Uout));

  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject*
ctimeIntegrationComputeElementNeighborShapeGradients(PyObject* self,
					      PyObject* args)
{
  PyObject *elementBoundariesArray, *elementNeighborsArray, *elementBarycentersArray,
    *elementBoundaryBarycentersArray, *elementNeighborShapeGradients;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &elementBoundariesArray,
                       &elementNeighborsArray,
                       &elementBarycentersArray,
                       &elementBoundaryBarycentersArray,
                       &elementNeighborShapeGradients))
    return NULL;

  computeElementNeighborShapeGradients(SHAPE(elementNeighborsArray)[0],
				       SHAPE(elementNeighborsArray)[1],
				       SHAPE(elementNeighborsArray)[1]-1,/*nSpace from simplex size*/
				       IDATA(elementBoundariesArray),
				       IDATA(elementNeighborsArray),
				       DDATA(elementBarycentersArray),
				       DDATA(elementBoundaryBarycentersArray),
				       DDATA(elementNeighborShapeGradients));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
ctimeIntegrationComputeCockburnDGlimiterArrays2d(PyObject* self,
					  PyObject* args)
{
  PyObject *elementBoundariesArray, *elementNeighborsArray, *elementBarycentersArray,
    *elementBoundaryBarycentersArray, *elementNeighborShapeGradients,
    *alphas, *alphaNeighbors;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &elementBoundariesArray,
                       &elementNeighborsArray,
                       &elementBarycentersArray,
                       &elementBoundaryBarycentersArray,
                       &elementNeighborShapeGradients,
		       &alphas,
		       &alphaNeighbors))
    return NULL;

  computeCockburnDGlimiterArrays2d(SHAPE(elementNeighborsArray)[0],
				   SHAPE(elementNeighborsArray)[1],
				   SHAPE(elementNeighborsArray)[1]-1,/*nSpace from simplex size*/
				   IDATA(elementBoundariesArray),
				   IDATA(elementNeighborsArray),
				   DDATA(elementBarycentersArray),
				   DDATA(elementBoundaryBarycentersArray),
				   DDATA(elementNeighborShapeGradients),
				   DDATA(alphas),
				   IDATA(alphaNeighbors));

  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
ctimeIntegrationApplyCockburnDGlimiterP1Lagrange2d(PyObject* self,
					    PyObject* args)
{
  PyObject *elementNeighborsArray, *l2g, *tag,
    *alphas, *alphaNeighbors, *Uin,*Uout;
  double nu,Mh2;
  if(!PyArg_ParseTuple(args,"ddOOOOOOO",
		       &nu,
		       &Mh2,
                       &elementNeighborsArray,
		       &l2g,
		       &tag,
		       &alphas,
		       &alphaNeighbors,
		       &Uin,
		       &Uout))
    return NULL;

  applyCockburnDGlimiterP1Lagrange2d(nu,Mh2,
				     SHAPE(elementNeighborsArray)[0],
				     SHAPE(elementNeighborsArray)[1],
				     SHAPE(elementNeighborsArray)[1]-1,/*nSpace from simplex size*/
				     SHAPE(l2g)[1],
				     IDATA(elementNeighborsArray),
				     IDATA(l2g),
				     IDATA(tag),
				     DDATA(alphas),
				     IDATA(alphaNeighbors),
				     DDATA(Uin),
				     DDATA(Uout));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
ctimeIntegrationApplyDurlofskyDGlimiterP1Lagrange2d(PyObject* self,
						    PyObject* args)
{
  PyObject *elementNeighborsArray, *elementBoundariesArray, 
    *elementNodesArray, *nodeArray,
    *elementBarycentersArray, *elementBoundaryBarycentersArray,  *elementNeighborShapeGradients,
    *l2g, *grad_v, *elementAverages, *tag, *Uin, *Uout;
  int killExtrema,allowMinWithUndershoot;
  if(!PyArg_ParseTuple(args,"iiOOOOOOOOOOOOO",
		       &killExtrema,
		       &allowMinWithUndershoot,
                       &elementNeighborsArray,
                       &elementBoundariesArray,
                       &elementNodesArray,
                       &nodeArray,
                       &elementBarycentersArray,
                       &elementBoundaryBarycentersArray,
                       &elementNeighborShapeGradients,
		       &l2g,
		       &grad_v,
		       &elementAverages,
		       &tag,
		       &Uin,
		       &Uout))
    return NULL;

  applyDurlofskyDGlimiterP1Lagrange2d(killExtrema,
				      allowMinWithUndershoot,
				      SHAPE(elementNeighborsArray)[0],
				      SHAPE(elementNeighborsArray)[1],
				      SHAPE(elementNodesArray)[1],
				      SHAPE(grad_v)[3],
				      SHAPE(l2g)[1],
				      IDATA(elementNeighborsArray),
				      IDATA(elementBoundariesArray),
				      IDATA(elementNodesArray),
				      DDATA(nodeArray),
				      DDATA(elementBarycentersArray),
				      DDATA(elementBoundaryBarycentersArray),
				      DDATA(elementNeighborShapeGradients),
				      IDATA(l2g),
				      DDATA(grad_v),
				      DDATA(elementAverages),
				      IDATA(tag),
				      DDATA(Uin),
				      DDATA(Uout));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
ctimeIntegrationApplyDurlofskyDGlimiterP1Lagrange2d_withVacuumTol(PyObject* self,
								  PyObject* args)
{
  PyObject *elementNeighborsArray, *elementBoundariesArray, 
    *elementNodesArray, *nodeArray,
    *elementBarycentersArray, *elementBoundaryBarycentersArray,  *elementNeighborShapeGradients,
    *l2g, *grad_v, *elementAverages, *tag, *Uin, *Uout;
  int killExtrema,allowMinWithUndershoot;
  int enforcePositivity;
  double vacuumTol;
  if(!PyArg_ParseTuple(args,"iiidOOOOOOOOOOOOO",
		       &killExtrema,
		       &allowMinWithUndershoot,
		       &enforcePositivity,
		       &vacuumTol,
                       &elementNeighborsArray,
                       &elementBoundariesArray,
                       &elementNodesArray,
                       &nodeArray,
                       &elementBarycentersArray,
                       &elementBoundaryBarycentersArray,
                       &elementNeighborShapeGradients,
		       &l2g,
		       &grad_v,
		       &elementAverages,
		       &tag,
		       &Uin,
		       &Uout))
    return NULL;

  applyDurlofskyDGlimiterP1Lagrange2d_withVacuumTol(killExtrema,
						    allowMinWithUndershoot,
						    enforcePositivity,
						    vacuumTol,
						    SHAPE(elementNeighborsArray)[0],
						    SHAPE(elementNeighborsArray)[1],
						    SHAPE(elementNodesArray)[1],
						    SHAPE(grad_v)[3],
						    SHAPE(l2g)[1],
						    IDATA(elementNeighborsArray),
						    IDATA(elementBoundariesArray),
						    IDATA(elementNodesArray),
						    DDATA(nodeArray),
						    DDATA(elementBarycentersArray),
						    DDATA(elementBoundaryBarycentersArray),
						    DDATA(elementNeighborShapeGradients),
						    IDATA(l2g),
						    DDATA(grad_v),
						    DDATA(elementAverages),
						    IDATA(tag),
						    DDATA(Uin),
						    DDATA(Uout));
  
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
ctimeIntegrationApplyDurlofskyDGlimiterP1Lagrange3d(PyObject* self,
					     PyObject* args)
{
  PyObject *elementNeighborsArray, *elementBoundariesArray, 
    *elementNodesArray, *nodeArray,
    *elementBarycentersArray, *elementBoundaryBarycentersArray,  *elementNeighborShapeGradients,
    *l2g, *grad_v, *elementAverages, *tag, *Uin, *Uout;
  int killExtrema,allowMinWithUndershoot;
  if(!PyArg_ParseTuple(args,"iiOOOOOOOOOOOOO",
		       &killExtrema,
		       &allowMinWithUndershoot,
                       &elementNeighborsArray,
                       &elementBoundariesArray,
                       &elementNodesArray,
                       &nodeArray,
                       &elementBarycentersArray,
                       &elementBoundaryBarycentersArray,
                       &elementNeighborShapeGradients,
		       &l2g,
		       &grad_v,
		       &elementAverages,
		       &tag,
		       &Uin,
		       &Uout))
    return NULL;

  applyDurlofskyDGlimiterP1Lagrange3d(killExtrema,
				      allowMinWithUndershoot,
				      SHAPE(elementNeighborsArray)[0],
				      SHAPE(elementNeighborsArray)[1],
				      SHAPE(elementNodesArray)[1],
				      SHAPE(grad_v)[3],
				      SHAPE(l2g)[1],
				      IDATA(elementNeighborsArray),
				      IDATA(elementBoundariesArray),
				      IDATA(elementNodesArray),
				      DDATA(nodeArray),
				      DDATA(elementBarycentersArray),
				      DDATA(elementBoundaryBarycentersArray),
				      DDATA(elementNeighborShapeGradients),
				      IDATA(l2g),
				      DDATA(grad_v),
				      DDATA(elementAverages),
				      IDATA(tag),
				      DDATA(Uin),
				      DDATA(Uout));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef ctimeIntegrationMethods[] = {
  { "psiTCtteDT",
    ctimeIntegrationPsiTCtteDT,
    METH_VARARGS,
    "Pseudo-transient continuation TTE time step selection strategy"},
  {"applyDGlimitingP1Lagrange1d",
   ctimeIntegrationApplyDGlimitingP1Lagrange1d,
   METH_VARARGS,
   "simple 1d DG limiting assuming P1 Lagrange local rep"},
  {"applyDGlimitingP1Lagrange1d_withVacuumTol",
   ctimeIntegrationApplyDGlimitingP1Lagrange1d_withVacuumTol,
   METH_VARARGS,
   "simple 1d DG limiting assuming P1 Lagrange local rep and enforcing more conservative limiting around u=0 state"},
  {"computeElementNeighborShapeGradients",
   ctimeIntegrationComputeElementNeighborShapeGradients,
   METH_VARARGS,
   "compute local shape functions using neighboring element barycenters"},
  {"computeCockburnDGlimiterArrays2d",
   ctimeIntegrationComputeCockburnDGlimiterArrays2d,
   METH_VARARGS,
   "compute alpha coefficients for dg limiter from Cockburn's notes"},
  {"applyCockburnDGlimiterP1Lagrange2d",
   ctimeIntegrationApplyCockburnDGlimiterP1Lagrange2d,
   METH_VARARGS,
   "apply 2d dg limiter from Cockburn's notes"},
  {"applyDurlofskyDGlimiterP1Lagrange2d",
   ctimeIntegrationApplyDurlofskyDGlimiterP1Lagrange2d,
   METH_VARARGS,
   "apply 2d dg limiter based on Durlofsky, Enquist and Osher JCP 92"},
  {"applyDurlofskyDGlimiterP1Lagrange2d_withVacuumTol",
   ctimeIntegrationApplyDurlofskyDGlimiterP1Lagrange2d_withVacuumTol,
   METH_VARARGS,
   "apply 2d dg limiter based on Durlofsky, Enquist and Osher JCP 92 for Sw equations"},
  {"applyDurlofskyDGlimiterP1Lagrange3d",
   ctimeIntegrationApplyDurlofskyDGlimiterP1Lagrange3d,
   METH_VARARGS,
   "apply 3d dg limiter based on Durlofsky, Enquist and Osher JCP 92"},
  { NULL,NULL,0,NULL}
};


PyMODINIT_FUNC initctimeIntegration(void)
{
  PyObject *m,*d;
  m = Py_InitModule("ctimeIntegration", ctimeIntegrationMethods);
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
