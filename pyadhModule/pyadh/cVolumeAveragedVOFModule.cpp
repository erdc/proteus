#include "Python.h"
#include "numpy/arrayobject.h"
#include "VolumeAveragedVOF.h"
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

static PyObject* cVolumeAveragedVOF_calculateResidual(PyObject* self,
					   PyObject* args)
{
  int nElements_global,
    offset_u,stride_u;
  int lag_shockCapturing;
  double alphaBDF,
    shockCapturingDiffusion,eps;
  PyObject *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *velocity,
    *porosity,
    *q_m,
    *q_u,
    *q_m_betaBDF,
    *cfl,
    *numDiff_u, 
    *numDiff_u_last, 
    *elementResidual_u, 
    *globalResidual;
  int nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *u_trial_ext,
    *u_grad_trial_ext,
    *velocity_ext,
    *porosity_ext,
    *n_ext,
    *isDOFBoundary_u,
    *bc_u_ext,
    *isFluxBoundary_u,
    *bc_flux_u_ext,
    *u_test_dS_ext,
    *phi_ext,
    *u_ext,
    *flux_ext;
  if (!PyArg_ParseTuple(args,
                        "dididOOOOOOOOOOOOOOOOiiOiOOOOOOOOOOOOOOOO",
			&eps,
                        &nElements_global,
			&alphaBDF,
			&lag_shockCapturing,
			&shockCapturingDiffusion,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&velocity,
			&porosity,
			&q_m,
			&q_u,
			&q_m_betaBDF,
			&cfl,
			&numDiff_u,
			&numDiff_u_last,
			&elementResidual_u, 
			&offset_u,
                        &stride_u,
			&globalResidual,
			&nExteriorElementBoundaries_global,
                        &exteriorElementBoundariesArray,
                        &elementBoundaryElementsArray,
                        &elementBoundaryLocalElementBoundariesArray,
                        &u_trial_ext,
                        &u_grad_trial_ext,
                        &velocity_ext,
			&porosity_ext,
                        &n_ext,
                        &isDOFBoundary_u,
                        &bc_u_ext,
                        &isFluxBoundary_u,
                        &bc_flux_u_ext,
                        &u_test_dS_ext,
			&phi_ext,
			&u_ext,
			&flux_ext))
    return NULL;
  
  //calculateResidual_VolumeAveragedVOF(nElements_global,
  VOL_VOF_RES(eps,
	      nElements_global,
	      alphaBDF,
	      lag_shockCapturing,
	      shockCapturingDiffusion,
	      IDATA(u_l2g), 
	      DDATA(elementDiameter),
	      DDATA(u_dof), 
	      DDATA(u_trial), 
	      DDATA(u_grad_trial), 
	      DDATA(u_test_dV), 
	      DDATA(u_grad_test_dV), 
	      DDATA(velocity),
	      DDATA(porosity),
	      DDATA(q_m),
	      DDATA(q_u),
	      DDATA(q_m_betaBDF),
	      DDATA(cfl),
	      DDATA(numDiff_u),
	      DDATA(numDiff_u_last),
	      DDATA(elementResidual_u), 
	      offset_u,stride_u,
	      DDATA(globalResidual),
	      nExteriorElementBoundaries_global,
	      IDATA(exteriorElementBoundariesArray),
	      IDATA(elementBoundaryElementsArray),
	      IDATA(elementBoundaryLocalElementBoundariesArray),
	      DDATA(u_trial_ext),
	      DDATA(u_grad_trial_ext),
	      DDATA(velocity_ext),
	      DDATA(porosity_ext),
	      DDATA(n_ext),
	      IDATA(isDOFBoundary_u),
	      DDATA(bc_u_ext),
	      IDATA(isFluxBoundary_u),
	      DDATA(bc_flux_u_ext),
	      DDATA(u_test_dS_ext),
	      DDATA(phi_ext),
	      DDATA(u_ext),
	      DDATA(flux_ext));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cVolumeAveragedVOF_calculateJacobian(PyObject* self,
					 PyObject* args)
{
  int nElements_global;
  int lag_shockCapturing;
  double alphaBDF,
    shockCapturingDiffusion;
  PyObject  *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *velocity,
    *porosity,
    *q_m_betaBDF, 
    *cfl,
    *numDiff_u_last; 
  PyObject *csrRowIndeces_u_u,*csrColumnOffsets_u_u;
  PyObject* globalJacobian;
  int nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *u_trial_ext,
    *u_grad_trial_ext,
    *velocity_ext,
    *porosity_ext,
    *n_ext,
    *isDOFBoundary_u,
    *bc_u_ext,
    *isFluxBoundary_u,
    *bc_flux_u_ext,
    *u_test_dS_ext;
  PyObject *csrColumnOffsets_eb_u_u;
  if (!PyArg_ParseTuple(args,
                        "ididOOOOOOOOOOOOOOOiOOOOOOOOOOOOOO",
                        &nElements_global,
			&alphaBDF,
			&lag_shockCapturing,
			&shockCapturingDiffusion,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&velocity,
			&porosity,
			&q_m_betaBDF,
			&cfl,
			&numDiff_u_last,
			&csrRowIndeces_u_u,
			&csrColumnOffsets_u_u,
			&globalJacobian,
			&nExteriorElementBoundaries_global,
                        &exteriorElementBoundariesArray,
                        &elementBoundaryElementsArray,
                        &elementBoundaryLocalElementBoundariesArray,
			&u_trial_ext,
			&u_grad_trial_ext,
                        &velocity_ext,
			&porosity_ext,
                        &n_ext,
                        &isDOFBoundary_u,
                        &bc_u_ext,
                        &isFluxBoundary_u,
                        &bc_flux_u_ext,
                        &u_test_dS_ext,
			&csrColumnOffsets_eb_u_u))
    return NULL;
  //calculateJacobian_VolumeAveragedVOF(nElements_global,
  VOL_VOF_JAC(nElements_global,
	      alphaBDF,
	      lag_shockCapturing,
	      shockCapturingDiffusion,
	      IDATA(u_l2g), 
	      DDATA(elementDiameter),
	      DDATA(u_dof), 
	      DDATA(u_trial), 
	      DDATA(u_grad_trial), 
	      DDATA(u_test_dV), 
	      DDATA(u_grad_test_dV), 
	      DDATA(velocity),
	      DDATA(porosity),
	      DDATA(q_m_betaBDF),
	      DDATA(cfl),
	      DDATA(numDiff_u_last),
	      IDATA(csrRowIndeces_u_u),IDATA(csrColumnOffsets_u_u),
	      CSRVAL(globalJacobian),
	      nExteriorElementBoundaries_global,
	      IDATA(exteriorElementBoundariesArray),
	      IDATA(elementBoundaryElementsArray),
	      IDATA(elementBoundaryLocalElementBoundariesArray),
	      DDATA( u_trial_ext),
	      DDATA( u_grad_trial_ext),
	      DDATA(velocity_ext),
	      DDATA(porosity_ext),
	      DDATA(n_ext),
	      IDATA(isDOFBoundary_u),
	      DDATA(bc_u_ext),
	      IDATA(isFluxBoundary_u),
	      DDATA(bc_flux_u_ext),
	      DDATA(u_test_dS_ext),
	      IDATA(csrColumnOffsets_eb_u_u));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyMethodDef cVolumeAveragedVOFMethods[] = {
 { "calculateResidual",
    cVolumeAveragedVOF_calculateResidual,
   METH_VARARGS, 
   "Calculate the global residual for the non-conservative level set equation"},
 { "calculateJacobian",
    cVolumeAveragedVOF_calculateJacobian,
   METH_VARARGS, 
   "Calculate the global Jacobian for the non-conservative level set equation"},
 { NULL,NULL,0,NULL}
};

extern "C"
{
PyMODINIT_FUNC initcVolumeAveragedVOF(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cVolumeAveragedVOF", cVolumeAveragedVOFMethods);
  d = PyModule_GetDict(m);
  import_array();
}
}//extern "C"
