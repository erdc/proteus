#include "Python.h"
#include "numpy/arrayobject.h"
#include "testStuffImpl.h"
/** \file testStuffImplModule.c
    \defgroup testStuffImplModule testStuffImpleModule
    \brief Python interface to testStuffImpl
    \@{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define NSF(p) ((NodeStarFactor *) p)

static PyObject*
testStuffAdvanceStageP1_C0_GLS_lump(PyObject* self,
				    PyObject* args)
{
  double dt;
  PyObject *elementDiameters,*n,*absDetJ,*sqrtDetG,*elementQuadratureWeights,
    *phiIn,*H,*dH,*r,*phiOut,*wOut;
  PyObject *l2g;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOdOO",
                       &elementDiameters,
		       &n,
		       &absDetJ,
		       &sqrtDetG,
		       &elementQuadratureWeights,
		       &l2g,
		       &phiIn,
		       &H,
		       &dH,
		       &r,
		       &dt,
		       &phiOut,
		       &wOut))
    return NULL;
  
  advanceStageP1_C0_GLS_lump(SHAPE(n)[0],
			     SHAPE(n)[1],
			     SHAPE(l2g)[1],
			     SHAPE(n)[3],
			     SHAPE(H)[1],
			     SHAPE(n)[2],
			     DDATA(elementDiameters),
			     DDATA(n),
			     DDATA(absDetJ),
			     DDATA(sqrtDetG),
			     DDATA(elementQuadratureWeights),
			     IDATA(l2g),
			     DDATA(phiIn),
			     DDATA(H),
			     DDATA(dH),
			     DDATA(r),
			     dt,
			     DDATA(phiOut),
			     DDATA(wOut));

  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject*
testStuffAdvanceStageP1_C0_GLS_lump_noSource(PyObject* self,
					     PyObject* args)
{
  double dt;
  PyObject *elementDiameters,*n,*absDetJ,*sqrtDetG,*elementQuadratureWeights,
    *phiIn,*H,*dH,*phiOut,*wOut;
  PyObject *l2g;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOdOO",
                       &elementDiameters,
		       &n,
		       &absDetJ,
		       &sqrtDetG,
		       &elementQuadratureWeights,
		       &l2g,
		       &phiIn,
		       &H,
		       &dH,
		       &dt,
		       &phiOut,
		       &wOut))
    return NULL;
  
  advanceStageP1_C0_GLS_lump_noSource(SHAPE(n)[0],
				      SHAPE(n)[1],
				      SHAPE(l2g)[1],
				      SHAPE(n)[3],
				      SHAPE(H)[1],
				      SHAPE(n)[2],
				      DDATA(elementDiameters),
				      DDATA(n),
				      DDATA(absDetJ),
				      DDATA(sqrtDetG),
				      DDATA(elementQuadratureWeights),
				      IDATA(l2g),
				      DDATA(phiIn),
				      DDATA(H),
				      DDATA(dH),
				      dt,
				      DDATA(phiOut),
				      DDATA(wOut));

  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject*
testStuffAdvanceStageP1_C0_SGS_lump_noSource(PyObject* self,
					     PyObject* args)
{
  double dt;
  PyObject *elementDiameters,*n,*absDetJ,*sqrtDetG,*elementQuadratureWeights,
    *phiIn,*H,*dH,*phiOut,*wOut;
  PyObject *l2g;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOdOO",
                       &elementDiameters,
		       &n,
		       &absDetJ,
		       &sqrtDetG,
		       &elementQuadratureWeights,
		       &l2g,
		       &phiIn,
		       &H,
		       &dH,
		       &dt,
		       &phiOut,
		       &wOut))
    return NULL;
  
  advanceStageP1_C0_SGS_lump_noSource(SHAPE(n)[0],
				      SHAPE(n)[1],
				      SHAPE(l2g)[1],
				      SHAPE(n)[3],
				      SHAPE(H)[1],
				      SHAPE(n)[2],
				      DDATA(elementDiameters),
				      DDATA(n),
				      DDATA(absDetJ),
				      DDATA(sqrtDetG),
				      DDATA(elementQuadratureWeights),
				      IDATA(l2g),
				      DDATA(phiIn),
				      DDATA(H),
				      DDATA(dH),
				      dt,
				      DDATA(phiOut),
				      DDATA(wOut));

  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject*
testStuffAdvanceStageP1_C0_SUPG_lump(PyObject* self,
				     PyObject* args)
{
  double dt;
  PyObject *elementDiameters,*n,*absDetJ,*sqrtDetG,*elementQuadratureWeights,
    *phiIn,*H,*dH,*r,*phiOut,*wOut;
  PyObject *l2g;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOdOO",
                       &elementDiameters,
		       &n,
		       &absDetJ,
		       &sqrtDetG,
		       &elementQuadratureWeights,
		       &l2g,
		       &phiIn,
		       &H,
		       &dH,
		       &r,
		       &dt,
		       &phiOut,
		       &wOut))
    return NULL;
  
  advanceStageP1_C0_SUPG_lump(SHAPE(n)[0],
			      SHAPE(n)[1],
			      SHAPE(l2g)[1],
			      SHAPE(n)[3],
			      SHAPE(H)[1],
			      SHAPE(n)[2],
			      DDATA(elementDiameters),
			      DDATA(n),
			      DDATA(absDetJ),
			      DDATA(sqrtDetG),
			      DDATA(elementQuadratureWeights),
			      IDATA(l2g),
			      DDATA(phiIn),
			      DDATA(H),
			      DDATA(dH),
			      DDATA(r),
			      dt,
			      DDATA(phiOut),
			      DDATA(wOut));

  Py_INCREF(Py_None);
  return Py_None;

}

static PyObject*
testStuffAdvanceStageRedistanceP1_C0_SUPG_lump(PyObject* self,
					       PyObject* args)
{
  double dt,eps;
  PyObject *elementDiameters,*n,*absDetJ,*sqrtDetG,*elementQuadratureWeights,
    *phiIn,*phiEvalIn,*phiOut,*wOut;
  PyObject *l2g;

  if(!PyArg_ParseTuple(args,"OOOOOOOOddOO",
                       &elementDiameters,
		       &n,
		       &absDetJ,
		       &sqrtDetG,
		       &elementQuadratureWeights,
		       &l2g,
		       &phiIn,
		       &phiEvalIn,
		       &dt,
		       &eps,
		       &phiOut,
		       &wOut))
    return NULL;
  
  advanceStageRedistanceP1_C0_SUPG_lump(SHAPE(n)[0],
					SHAPE(n)[1],
					SHAPE(l2g)[1],
					SHAPE(n)[3],
					SHAPE(elementQuadratureWeights)[0],
					SHAPE(n)[2],
					DDATA(elementDiameters),
					DDATA(n),
					DDATA(absDetJ),
					DDATA(sqrtDetG),
					DDATA(elementQuadratureWeights),
					IDATA(l2g),
					DDATA(phiIn),
					DDATA(phiEvalIn),
					dt,
					eps,
					DDATA(phiOut),
					DDATA(wOut));

  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject*
testStuffAdvanceStageRedistanceWeakDirP1_C0_SUPG_lump(PyObject* self,
						      PyObject* args)
{
  double dt,eps;
  PyObject *elementDiameters,*n,*absDetJ,*sqrtDetG,*elementQuadratureWeights,
    *phiIn,*phiEvalIn,*phiOut,*wOut;
  PyObject *l2g,*weakDirichletFlag,
    *exteriorElementBoundariesArray,*elementBoundaryElementsArray;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOddOOO",
                       &elementDiameters,
		       &n,
		       &absDetJ,
		       &sqrtDetG,
		       &elementQuadratureWeights,
		       &l2g,
		       &exteriorElementBoundariesArray,
		       &elementBoundaryElementsArray,
		       &phiIn,
		       &phiEvalIn,
		       &dt,
		       &eps,
		       &weakDirichletFlag,
		       &phiOut,
		       &wOut))
    return NULL;
  
  advanceStageRedistanceWeakDirP1_C0_SUPG_lump(SHAPE(n)[0],
					       SHAPE(n)[1],
					       SHAPE(l2g)[1],
					       SHAPE(n)[3],
					       SHAPE(elementQuadratureWeights)[0],
					       SHAPE(n)[2],
					       SHAPE(exteriorElementBoundariesArray)[0],
					       DDATA(elementDiameters),
					       DDATA(n),
					       DDATA(absDetJ),
					       DDATA(sqrtDetG),
					       DDATA(elementQuadratureWeights),
					       IDATA(l2g),
					       IDATA(exteriorElementBoundariesArray),
					       IDATA(elementBoundaryElementsArray),
					       DDATA(phiIn),
					       DDATA(phiEvalIn),
					       dt,
					       eps,
					       IDATA(weakDirichletFlag),
					       DDATA(phiOut),
					       DDATA(wOut));

  Py_INCREF(Py_None);
  return Py_None;

}


static PyMethodDef testStuffImplMethods[] = {
  { "advanceStageP1_C0_GLS_lump",
    testStuffAdvanceStageP1_C0_GLS_lump,
    METH_VARARGS,
    "Barth and Sethian stage update for P^1 C^0 affine elements with GLS type stab."},
  { "advanceStageP1_C0_GLS_lump_noSource",
    testStuffAdvanceStageP1_C0_GLS_lump_noSource,
    METH_VARARGS,
    "stage update for P^1 C^0 affine elements with GLS type stab and no source term"},
  { "advanceStageP1_C0_SGS_lump_noSource",
    testStuffAdvanceStageP1_C0_SGS_lump_noSource,
    METH_VARARGS,
    "stage update for P^1 C^0 affine elements with SGS type stab and no source term (same as GLS more or less)"},
  { "advanceStageP1_C0_SUPG_lump",
    testStuffAdvanceStageP1_C0_SUPG_lump,
    METH_VARARGS,
    "stage update for P^1 C^0 affine elements with SUPG type stab and source term "},
  { "advanceStageRedistanceP1_C0_SUPG_lump",
    testStuffAdvanceStageRedistanceP1_C0_SUPG_lump,
    METH_VARARGS,
    "redistancing stage update for P^1 C^0 affine elements with SUPG type stab and source term "},
  { "advanceStageRedistanceWeakDirP1_C0_SUPG_lump",
    testStuffAdvanceStageRedistanceWeakDirP1_C0_SUPG_lump,
    METH_VARARGS,
    "redistancing stage update for P^1 C^0 affine elements with SUPG type stab and source term \
     and attempt to freeze elements around zero level set "},
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC inittestStuffImpl(void)
{
  PyObject *m,*d;
  m = Py_InitModule3("testStuffImpl", testStuffImplMethods,"testStuff module");
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
