#include "Python.h"
#include "numpy/arrayobject.h"
#include "MCorrV2.h"
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

static PyObject* cMCorrV2_calculateResidual(PyObject* self,
					   PyObject* args)
{
  int nElements_global,
    offset_u,stride_u;
  double epsHeaviside,epsDirac,epsDiffusion;
  PyObject //testing
    *mesh_trial_ref,
    *mesh_grad_trial_ref,
    *mesh_dof,
    *mesh_l2g,
    *dV_ref,
    *u_trial_ref,
    *u_grad_trial_ref,
    *u_test_ref,
    *u_grad_test_ref,
    *mesh_trial_trace_ref,
    *mesh_grad_trial_trace_ref,
    *dS_ref,
    *u_trial_trace_ref,
    *u_grad_trial_trace_ref,
    *u_test_trace_ref,
    *u_grad_test_trace_ref,
    *normal_ref,
    *boundaryJac_ref,
    //testing
    *u_l2g, 
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
                        "OOOOOOOOOOOOOOOOOOidddOOOOOOOOOOOOiiO",
                        //testing
			&mesh_trial_ref,
			&mesh_grad_trial_ref,
			&mesh_dof,
			&mesh_l2g,
			&dV_ref,
			&u_trial_ref,
			&u_grad_trial_ref,
			&u_test_ref,
			&u_grad_test_ref,
			&mesh_trial_trace_ref,
			&mesh_grad_trial_trace_ref,
			&dS_ref,
			&u_trial_trace_ref,
			&u_grad_trial_trace_ref,
			&u_test_trace_ref,
			&u_grad_test_trace_ref,
			&normal_ref,
			&boundaryJac_ref,
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
  
  //calculateResidual_MCorrV2(nElements_global,
  MCORRV2_RES(//testing mesh replacement
			   DDATA(mesh_trial_ref),
			   DDATA(mesh_grad_trial_ref),
			   DDATA(mesh_dof),
			   IDATA(mesh_l2g),
			   DDATA(dV_ref),
			   DDATA(u_trial_ref),
			   DDATA(u_grad_trial_ref),
			   DDATA(u_test_ref),
			   DDATA(u_grad_test_ref),
			   DDATA(mesh_trial_trace_ref),
			   DDATA(mesh_grad_trial_trace_ref),
			   DDATA(dS_ref),
			   DDATA(u_trial_trace_ref),
			   DDATA(u_grad_trial_trace_ref),
			   DDATA(u_test_trace_ref),
			   DDATA(u_grad_test_trace_ref),
			   DDATA(normal_ref),
			   DDATA(boundaryJac_ref),
			   //end testing meshreplacement
			   nElements_global,
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

static PyObject* cMCorrV2_calculateJacobian(PyObject* self,
					 PyObject* args)
{
  int nElements_global;
  double epsHeaviside,epsDirac,epsDiffusion;
  PyObject  //testing
    *mesh_trial_ref,
    *mesh_grad_trial_ref,
    *mesh_dof,
    *mesh_l2g,
    *dV_ref,
    * u_trial_ref,
    * u_grad_trial_ref,
    * u_test_ref,
    * u_grad_test_ref,
    *mesh_trial_trace_ref,
    *mesh_grad_trial_trace_ref,
    *dS_ref,
    *u_trial_trace_ref,
    *u_grad_trial_trace_ref,
    *u_test_trace_ref,
    *u_grad_test_trace_ref,
    *normal_ref,
    *boundaryJac_ref,
    //testing
    *u_l2g, 
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
                        "OOOOOOOOOOOOOOOOOOidddOOOOOOOOOOOO",
                        //testing
			&mesh_trial_ref,
			&mesh_grad_trial_ref,
			&mesh_dof,
			&mesh_l2g,
			&dV_ref,
			&u_trial_ref,
			&u_grad_trial_ref,
			&u_test_ref,
			&u_grad_test_ref,
			&mesh_trial_trace_ref,
			&mesh_grad_trial_trace_ref,
			&dS_ref,
			&u_trial_trace_ref,
			&u_grad_trial_trace_ref,
			&u_test_trace_ref,
			&u_grad_test_trace_ref,
			&normal_ref,
			&boundaryJac_ref,
			//testing
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
  //calculateJacobian_MCorrV2(nElements_global,
  MCORRV2_JAC(//testing mesh replacement
  			     DDATA(mesh_trial_ref),
  			     DDATA(mesh_grad_trial_ref),
  			     DDATA(mesh_dof),
  			     IDATA(mesh_l2g),
  			     DDATA(dV_ref),
  			     DDATA(u_trial_ref),
  			     DDATA(u_grad_trial_ref),
  			     DDATA(u_test_ref),
  			     DDATA(u_grad_test_ref),
  			     DDATA(mesh_trial_trace_ref),
  			     DDATA(mesh_grad_trial_trace_ref),
  			     DDATA(dS_ref),
  			     DDATA(u_trial_trace_ref),
  			     DDATA(u_grad_trial_trace_ref),
  			     DDATA(u_test_trace_ref),
  			     DDATA(u_grad_test_trace_ref),
  			     DDATA(normal_ref),
  			     DDATA(boundaryJac_ref),
  			     //end testing meshreplacement
  			     nElements_global,
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

static PyMethodDef cMCorrV2Methods[] = {
 { "calculateResidual",
    cMCorrV2_calculateResidual,
   METH_VARARGS, 
   "Calculate the global residual for the non-conservative level set equation"},
 { "calculateJacobian",
    cMCorrV2_calculateJacobian,
   METH_VARARGS, 
   "Calculate the global Jacobian for the non-conservative level set equation"},
 { NULL,NULL,0,NULL}
};

extern "C"
{
PyMODINIT_FUNC initcMCorrV2(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cMCorrV2", cMCorrV2Methods);
  d = PyModule_GetDict(m);
  import_array();
}
}//extern "C"
