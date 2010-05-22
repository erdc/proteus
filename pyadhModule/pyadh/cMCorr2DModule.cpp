#include "Python.h"
#include "numpy/arrayobject.h"
#include "MCorr2D.h"
#include "superluWrappersModule.h"

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

static PyObject* cMCorr2D_calculateResidual(PyObject* self,
					   PyObject* args)
{
  int nElements_global,
    offset_u,stride_u;
  double epsHeaviside,epsDirac,epsDiffusion;
  PyObject *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *q_phi,
    *q_H,
    *q_u,
    *q_r,
    *elementResidual_u, 
    *globalResidual;
  if (!PyArg_ParseTuple(args,
                        "idddOOOOOOOOOOOOiiO",
                        &nElements_global,
			&epsHeaviside,
			&epsDirac,
			&epsDiffusion,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&q_phi,
			&q_H,
			&q_u,
			&q_r,
			&elementResidual_u, 
			&offset_u,&stride_u,
			&globalResidual))
    return NULL;
  
  //calculateResidual_MCorr2D(nElements_global,
  MCORR_RES(nElements_global,
		    epsHeaviside,
		    epsDirac,
		    epsDiffusion,
		    IDATA(u_l2g), 
		    DDATA(elementDiameter),
		    DDATA(u_dof), 
		    DDATA(u_trial), 
		    DDATA(u_grad_trial), 
		    DDATA(u_test_dV), 
		    DDATA(u_grad_test_dV), 
		    DDATA(q_phi),
		    DDATA(q_H),
		    DDATA(q_u),
		    DDATA(q_r),
		    offset_u,stride_u,
		    DDATA(globalResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cMCorr2D_calculateJacobian(PyObject* self,
					 PyObject* args)
{
  int nElements_global;
  double epsHeaviside,epsDirac,epsDiffusion;
  PyObject  *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *q_phi,
    *q_H,
    *csrRowIndeces_u_u,*csrColumnOffsets_u_u,
    *globalJacobian;
  if (!PyArg_ParseTuple(args,
                        "idddOOOOOOOOOOOO",
                        &nElements_global,
			&epsHeaviside,
			&epsDirac,
			&epsDiffusion,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&q_phi,
			&q_H,
			&csrRowIndeces_u_u,
			&csrColumnOffsets_u_u,
			&globalJacobian))
    return NULL;
  //calculateJacobian_MCorr2D(nElements_global,
  MCORR_JAC(nElements_global,
		    epsHeaviside,
		    epsDirac,
		    epsDiffusion,
		    IDATA(u_l2g), 
		    DDATA(elementDiameter),
		    DDATA(u_dof), 
		    DDATA(u_trial), 
		    DDATA(u_grad_trial), 
		    DDATA(u_test_dV), 
		    DDATA(u_grad_test_dV), 
		    DDATA(q_phi),
		    DDATA(q_H),
		    IDATA(csrRowIndeces_u_u),IDATA(csrColumnOffsets_u_u),
		    CSRVAL(globalJacobian));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef cMCorr2DMethods[] = {
 { "calculateResidual",
    cMCorr2D_calculateResidual,
   METH_VARARGS, 
   "Calculate the global residual for the non-conservative level set equation"},
 { "calculateJacobian",
    cMCorr2D_calculateJacobian,
   METH_VARARGS, 
   "Calculate the global Jacobian for the non-conservative level set equation"},
 { NULL,NULL,0,NULL}
};

extern "C"
{
PyMODINIT_FUNC initcMCorr2D(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cMCorr2D", cMCorr2DMethods);
  d = PyModule_GetDict(m);
  import_array();
}
}//extern "C"
