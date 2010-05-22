#include "Python.h"
#include "numpy/arrayobject.h"
#include "RDLSQ.h"
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

static PyObject* cRDLSQ_calculateResidual(PyObject* self,
					   PyObject* args)
{
  int nElements_global,
    offset_u,stride_u;
  int useTimeIntegration,lag_shockCapturing,lag_stabilization,
    freezeLevelSet;
  double dt,eps,
    shockCapturingDiffusion;
  PyObject *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *phi_ls,
    *q_m,
    *q_u,
    *q_dH,
    *u_weak_internal_bc_dofs,
    *q_m_last,
    *q_dH_last,
    *cfl,
    *numDiff_u, 
    *numDiff_u_last, 
    *elementResidual_u,
    //for debugging
    *q_m_t,
    *q_r,
    *q_subgridError,
    *q_Lstar,
    *q_tau_last,
    //mwf end debugging 
    *weakDirichletConditionFlags,
    *globalResidual;
  int nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *u_trial_ext,
    *u_grad_trial_ext,
    *phi_ls_ext,
    *n_ext,
    *isDOFBoundary_u,
    *bc_u_ext,
    *u_test_dS_ext,
    *u_ext;
  if (!PyArg_ParseTuple(args,
                        "iddiiiidOOOOOOOOOOOOOOOOOOOOOOOOiiOiOOOOOOOOOOO",//added 4*O for debugging
                        &nElements_global,
			&dt,
			&eps,
			&freezeLevelSet,
			&useTimeIntegration,
			&lag_shockCapturing,
			&lag_stabilization,
			&shockCapturingDiffusion,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&phi_ls,
			&q_m,
			&q_u,
			&q_dH,
			&u_weak_internal_bc_dofs,
			&q_m_last,
			&q_dH_last,
			&cfl,
			&numDiff_u,
			&numDiff_u_last,
			&elementResidual_u,
			//mwf for debugging
			&q_m_t,
			&q_r,
			&q_subgridError,
			&q_Lstar,
			&q_tau_last,
			//mwf end debugging 
			&weakDirichletConditionFlags,
			&offset_u,&stride_u,
			&globalResidual,
			&nExteriorElementBoundaries_global,
                        &exteriorElementBoundariesArray,
                        &elementBoundaryElementsArray,
                        &elementBoundaryLocalElementBoundariesArray,
                        &u_trial_ext,
                        &u_grad_trial_ext,
                        &phi_ls_ext,
                        &n_ext,
                        &isDOFBoundary_u,
                        &bc_u_ext,
                        &u_test_dS_ext,
			&u_ext))
    return NULL;
  
//   calculateResidual_RDLSQ(nElements_global,
  RDLSQ_RES (nElements_global,
			 dt,
			 eps,
			 freezeLevelSet,
			 useTimeIntegration,
			 lag_shockCapturing,
			 lag_stabilization,
			 shockCapturingDiffusion,
			 IDATA(u_l2g), 
			 DDATA(elementDiameter),
			 DDATA(u_dof), 
			 DDATA(u_trial), 
			 DDATA(u_grad_trial), 
			 DDATA(u_test_dV), 
			 DDATA(u_grad_test_dV), 
			 DDATA(phi_ls),
			 DDATA(q_m),
			 DDATA(q_u),
			 DDATA(q_dH),
			 DDATA(u_weak_internal_bc_dofs),
			 DDATA(q_m_last),
			 DDATA(q_dH_last),
			 DDATA(cfl),
			 DDATA(numDiff_u),
			 DDATA(numDiff_u_last),
			 DDATA(elementResidual_u), 
			 //mwf for debugging
			 DDATA(q_m_t),
			 DDATA(q_r),
			 DDATA(q_subgridError),
			 DDATA(q_Lstar),
			 DDATA(q_tau_last),
			 //mwf end debugging 
			 IDATA(weakDirichletConditionFlags),
			 offset_u,stride_u,
			 DDATA(globalResidual),
			 nExteriorElementBoundaries_global,
			 IDATA(exteriorElementBoundariesArray),
			 IDATA(elementBoundaryElementsArray),
			 IDATA(elementBoundaryLocalElementBoundariesArray),
			 DDATA(u_trial_ext),
			 DDATA(u_grad_trial_ext),
			 DDATA(phi_ls_ext),
			 DDATA(n_ext),
			 IDATA(isDOFBoundary_u),
			 DDATA(bc_u_ext),
			 DDATA(u_test_dS_ext),
			 DDATA(u_ext));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cRDLSQ_calculateJacobian(PyObject* self,
					 PyObject* args)
{
  int nElements_global;
  int useTimeIntegration,lag_shockCapturing,lag_stabilization,
    freezeLevelSet;
  double dt,eps,
    shockCapturingDiffusion;
  PyObject  *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *phi_ls,
    *weakDirichletConditionFlags,
    *q_m_last,
    *q_dH_last,
    *cfl,
    *numDiff_u,
    *numDiff_u_last; 
  PyObject *csrRowIndeces_u_u,*csrColumnOffsets_u_u;
  PyObject* globalJacobian;
  int nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *u_trial_ext,
    *u_grad_trial_ext,
    *phi_ls_ext,
    *n_ext,
    *isDOFBoundary_u,
    *bc_u_ext,
    *u_test_dS_ext;
  PyObject *csrColumnOffsets_eb_u_u;
  if (!PyArg_ParseTuple(args,
                        "iddiiiidOOOOOOOOOOOOOOOOOiOOOOOOOOOOO",
                        &nElements_global,
			&dt,
			&eps,
			&freezeLevelSet,
			&useTimeIntegration,
			&lag_shockCapturing,
			&lag_stabilization,
			&shockCapturingDiffusion,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&phi_ls,
			&q_m_last,
			&q_dH_last,
			&cfl,
			&numDiff_u,
			&numDiff_u_last,
			&weakDirichletConditionFlags,
			&csrRowIndeces_u_u,
			&csrColumnOffsets_u_u,
			&globalJacobian,
			&nExteriorElementBoundaries_global,
                        &exteriorElementBoundariesArray,
                        &elementBoundaryElementsArray,
                        &elementBoundaryLocalElementBoundariesArray,
			&u_trial_ext,
			&u_grad_trial_ext,
                        &phi_ls_ext,
                        &n_ext,
                        &isDOFBoundary_u,
                        &bc_u_ext,
                        &u_test_dS_ext,
			&csrColumnOffsets_eb_u_u))
    return NULL;
  //  calculateJacobian_RDLSQ(nElements_global, 
  RDLSQ_JAC (nElements_global, 
			 dt,
			 eps,
			 freezeLevelSet,
			 useTimeIntegration,
			 lag_shockCapturing,
			 lag_stabilization,
			 shockCapturingDiffusion,
			 IDATA(u_l2g), 
			 DDATA(elementDiameter),
			 DDATA(u_dof), 
			 DDATA(u_trial), 
			 DDATA(u_grad_trial), 
			 DDATA(u_test_dV), 
			 DDATA(u_grad_test_dV), 
			 DDATA(phi_ls),
			 DDATA(q_m_last),
			 DDATA(q_dH_last),
			 DDATA(cfl),
			 DDATA(numDiff_u),
			 DDATA(numDiff_u_last),
			 IDATA(weakDirichletConditionFlags),
			 IDATA(csrRowIndeces_u_u),IDATA(csrColumnOffsets_u_u),
			 CSRVAL(globalJacobian),
			 nExteriorElementBoundaries_global,
			 IDATA(exteriorElementBoundariesArray),
			 IDATA(elementBoundaryElementsArray),
			 IDATA(elementBoundaryLocalElementBoundariesArray),
			 DDATA( u_trial_ext),
			 DDATA( u_grad_trial_ext),
			 DDATA(phi_ls_ext),
			 DDATA(n_ext),
			 IDATA(isDOFBoundary_u),
			 DDATA(bc_u_ext),
			 DDATA(u_test_dS_ext),
			 IDATA(csrColumnOffsets_eb_u_u));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyMethodDef cRDLSQMethods[] = {
 { "calculateResidual",
    cRDLSQ_calculateResidual,
   METH_VARARGS, 
   "Calculate the global residual for redistancing the level set equation"},
 { "calculateJacobian",
    cRDLSQ_calculateJacobian,
   METH_VARARGS, 
   "Calculate the global Jacobian for redistancing the level set equation"},
 { NULL,NULL,0,NULL}
};

extern "C"
{
PyMODINIT_FUNC initcRDLSQ(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cRDLSQ", cRDLSQMethods);
  d = PyModule_GetDict(m);
  import_array();
}
}//extern "C"
