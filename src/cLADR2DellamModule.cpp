#include "Python.h"
#include "numpy/arrayobject.h"
#include "LADR2Dellam.h"
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

static PyObject* cLADR2Dellam_calculateResidual(PyObject* self,
						PyObject* args)
{
  int nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nQuadraturePoints_element,
    nQuadraturePoints_elementBoundary,
    offset_u,stride_u;
  double theta,aL,aT,Dm,dtnp1;

  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *elementBoundaryOuterNormalsArray,
    *dV,
    *x_track,
    *t_track,
    *element_track,
    *flag_track,
    *x_dt,
    *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *velocity,
    *q_u,
    *q_m,
    *q_m_last,
    *cfl,
    *elementResidual_u, 
    *globalResidual,
    *sdInfo_u_rowptr,
    *sdInfo_u_colind;
  int nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *u_trial_ext,
    *u_grad_trial_ext,
    *velocity_ext,
    *n_ext,
    *isDOFBoundary_u,
    *bc_u_ext,
    *isFluxBoundary_u,
    *bc_flux_u_ext,
    *ebqe_outflow_flux_last,
    *u_test_dS_ext,
    *u_ext,
    *ebqe_outflow_flux;
  if (!PyArg_ParseTuple(args,
                        "dddddiiiiiOOOOOOOOOOOOOOOOOOOOOOOOOiiOiiOOOOOOOOOOOOOOO",
                        &theta,
			&aL,
			&aT,
			&Dm,
			&dtnp1,
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nQuadraturePoints_element,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&elementBoundaryOuterNormalsArray,
			&dV,
			&x_track,
			&t_track,
			&element_track,
			&flag_track,
			&x_dt,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&velocity,
			&q_u,
			&q_m,
			&q_m_last,
			&cfl,
			&elementResidual_u,
			&sdInfo_u_rowptr,
			&sdInfo_u_colind,
			&offset_u,
			&stride_u,
			&globalResidual,
			&nExteriorElementBoundaries_global,
			&nQuadraturePoints_elementBoundary,
                        &exteriorElementBoundariesArray,
                        &elementBoundaryElementsArray,
                        &elementBoundaryLocalElementBoundariesArray,
                        &u_trial_ext,
                        &u_grad_trial_ext,
                        &velocity_ext,
                        &n_ext,
                        &isDOFBoundary_u,
                        &bc_u_ext,
                        &isFluxBoundary_u,
                        &bc_flux_u_ext,
			&ebqe_outflow_flux_last,
                        &u_test_dS_ext,
			&u_ext,
			&ebqe_outflow_flux))
    return NULL;
  
  calculateResidual_LADR2Dellam(theta,
				aL,
				aT,
				Dm,
				dtnp1,
				nElements_global,
				nNodes_global,
				nNodes_element,
				nElementBoundaries_element,
				nQuadraturePoints_element,
				DDATA(nodeArray),
				IDATA(elementNodesArray),
				IDATA(elementNeighborsArray),
				DDATA(elementBoundaryOuterNormalsArray),
				DDATA(dV),
				DDATA(x_track),
				DDATA(t_track),
				IDATA(element_track),
				IDATA(flag_track),
				DDATA(x_dt),
				IDATA(u_l2g), 
				DDATA(elementDiameter),
				DDATA(u_dof), 
				DDATA(u_trial), 
				DDATA(u_grad_trial), 
				DDATA(u_test_dV), 
				DDATA(u_grad_test_dV), 
				DDATA(velocity),
				DDATA(q_u),
				DDATA(q_m),
				DDATA(q_m_last),
				DDATA(cfl),
				DDATA(elementResidual_u), 
				IDATA(sdInfo_u_rowptr),
				IDATA(sdInfo_u_colind),
				offset_u,stride_u,
				DDATA(globalResidual),
				nExteriorElementBoundaries_global,
				nQuadraturePoints_elementBoundary,
				IDATA(exteriorElementBoundariesArray),
				IDATA(elementBoundaryElementsArray),
				IDATA(elementBoundaryLocalElementBoundariesArray),
				DDATA(u_trial_ext),
				DDATA(u_grad_trial_ext),
				DDATA(velocity_ext),
				DDATA(n_ext),
				IDATA(isDOFBoundary_u),
				DDATA(bc_u_ext),
				IDATA(isFluxBoundary_u),
				DDATA(bc_flux_u_ext),
				DDATA(ebqe_outflow_flux_last),
				DDATA(u_test_dS_ext),
				DDATA(u_ext),
				DDATA(ebqe_outflow_flux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cLADR2Dellam_calculateJacobian(PyObject* self,
						PyObject* args)
{
  int nElements_global,
    nQuadraturePoints_element,nQuadraturePoints_elementBoundary;
  double theta,aL,aT,Dm,dtnp1;

  PyObject  *x_dt,
    *u_l2g, 
    *elementDiameter,
    *u_dof, 
    *u_trial, 
    *u_grad_trial, 
    *u_test_dV, 
    *u_grad_test_dV, 
    *velocity,
    *q_m_last, 
    *cfl,
    *sdInfo_u_rowptr, 
    *sdInfo_u_colind; 
  PyObject *csrRowIndeces_u_u,*csrColumnOffsets_u_u;
  PyObject* globalJacobian;
  int nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *u_trial_ext,
    *u_grad_trial_ext,
    *velocity_ext,
    *n_ext,
    *isDOFBoundary_u,
    *bc_u_ext,
    *isFluxBoundary_u,
    *bc_flux_u_ext,
    *u_test_dS_ext;
  PyObject *csrColumnOffsets_eb_u_u;
  if (!PyArg_ParseTuple(args,
                        "dddddiiOOOOOOOOOOOOOOOOiiOOOOOOOOOOOOO",
			&theta,
			&aL,
			&aT,
			&Dm,
			&dtnp1,
                        &nElements_global,
			&nQuadraturePoints_element,
			&x_dt,
			&u_l2g, 
			&elementDiameter,
			&u_dof, 
			&u_trial, 
			&u_grad_trial, 
			&u_test_dV, 
			&u_grad_test_dV, 
			&velocity,
			&q_m_last,
			&cfl,
			&sdInfo_u_rowptr,
			&sdInfo_u_colind,
			&csrRowIndeces_u_u,
			&csrColumnOffsets_u_u,
			&globalJacobian,
			&nExteriorElementBoundaries_global,
			&nQuadraturePoints_elementBoundary,
                        &exteriorElementBoundariesArray,
                        &elementBoundaryElementsArray,
                        &elementBoundaryLocalElementBoundariesArray,
			&u_trial_ext,
			&u_grad_trial_ext,
                        &velocity_ext,
                        &n_ext,
                        &isDOFBoundary_u,
                        &bc_u_ext,
                        &isFluxBoundary_u,
                        &bc_flux_u_ext,
                        &u_test_dS_ext,
			&csrColumnOffsets_eb_u_u))
    return NULL;
  calculateJacobian_LADR2Dellam(theta,
				aL,
				aT,
				Dm,
				dtnp1,
				nElements_global,
				nQuadraturePoints_element,
				DDATA(x_dt),
				IDATA(u_l2g), 
				DDATA(elementDiameter),
				DDATA(u_dof), 
				DDATA(u_trial), 
				DDATA(u_grad_trial), 
				DDATA(u_test_dV), 
				DDATA(u_grad_test_dV), 
				DDATA(velocity),
				DDATA(q_m_last),
				DDATA(cfl),
				IDATA(sdInfo_u_rowptr),
				IDATA(sdInfo_u_colind),
				IDATA(csrRowIndeces_u_u),IDATA(csrColumnOffsets_u_u),
				CSRVAL(globalJacobian),
				nExteriorElementBoundaries_global,
				nQuadraturePoints_elementBoundary,
				IDATA(exteriorElementBoundariesArray),
				IDATA(elementBoundaryElementsArray),
				IDATA(elementBoundaryLocalElementBoundariesArray),
				DDATA( u_trial_ext),
				DDATA( u_grad_trial_ext),
				DDATA(velocity_ext),
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
static PyObject* cLADR2Dellam_markInflowBoundaryPoints(PyObject* self,
						       PyObject* args)
{
  double tn,tnp1,t;

  int nExteriorElementBoundaries_global,
    nQuadraturePoints_elementBoundary;

  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *elementBoundaryOuterNormalsArray,
    *ebqe_x,
    *ebqe_n_ext,
    *ebqe_velocity_ext_last,
    *ebqe_velocity_ext,
    *isDOFBoundary_u,
    *isFluxBoundary_u,
    *element_track,
    *flag_track;

  if (!PyArg_ParseTuple(args,
                        "dddiiOOOOOOOOOOO",
			&tn,
			&tnp1,
			&t,
			&nExteriorElementBoundaries_global,
			&nQuadraturePoints_elementBoundary,
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
			&ebqe_x,
			&ebqe_n_ext,
			&ebqe_velocity_ext_last,
			&ebqe_velocity_ext,
			&isDOFBoundary_u,
			&isFluxBoundary_u,
			&element_track,
			&flag_track))
    return NULL;

  markInflowBoundaryPoints2d(tn,
			   tnp1,
			   t,
			   nExteriorElementBoundaries_global,
			   nQuadraturePoints_elementBoundary,
			   IDATA(exteriorElementBoundariesArray),
			   IDATA(elementBoundaryElementsArray),
			   IDATA(elementBoundaryLocalElementBoundariesArray),	
			     DDATA(ebqe_x),
			   DDATA(ebqe_n_ext),
			   DDATA(ebqe_velocity_ext_last),
			   DDATA(ebqe_velocity_ext),
			   IDATA(isDOFBoundary_u),
			   IDATA(isFluxBoundary_u),
			   IDATA(element_track),
			   IDATA(flag_track));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cLADR2Dellam_accumulateInflowFluxInGlobalResidual(PyObject* self,
								   PyObject* args)
{
  double tp,timeWeight;

  int nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    offset_u,stride_u,
    nExteriorElementBoundaries_global,
    nQuadraturePoints_elementBoundary;

  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray, 
    *elementBoundaryOuterNormalsArray,
    *dS,
    *x_track_ext,
    *t_track_ext,
    *element_track_ext,
    *flag_track_ext,
    *u_l2g,
    *u_dof,
    *elementResidual_u, 
    *globalResidual,
    *sdInfo_u_rowptr,
    *sdInfo_u_colind,
    *isFluxBoundary_u,
    *ebqe_bc_flux_u_ext;

  if (!PyArg_ParseTuple(args,
			"iiiiiiOOOOOOOddOOOOOOOOOOiiOOO",
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nExteriorElementBoundaries_global,
			&nQuadraturePoints_elementBoundary,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray, 
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
			&elementBoundaryOuterNormalsArray,
			&tp,
			&timeWeight,
			&dS,
			&x_track_ext,
			&t_track_ext,
			&element_track_ext,
			&flag_track_ext,
			&u_l2g, 
			&u_dof,
			&elementResidual_u, 
			&sdInfo_u_rowptr,
			&sdInfo_u_colind,
			&offset_u,
			&stride_u, 
			&isFluxBoundary_u,
			&ebqe_bc_flux_u_ext,
			&globalResidual))
    return NULL;

  accumulateInflowFluxInGlobalResidual2d(nElements_global,
				       nNodes_global,
				       nNodes_element,
				       nElementBoundaries_element,
				       nExteriorElementBoundaries_global,
				       nQuadraturePoints_elementBoundary,
				       DDATA(nodeArray),
				       IDATA(elementNodesArray),
				       IDATA(elementNeighborsArray), 
				       IDATA(exteriorElementBoundariesArray),
				       IDATA(elementBoundaryElementsArray),
				       IDATA(elementBoundaryLocalElementBoundariesArray),
				       DDATA(elementBoundaryOuterNormalsArray),
				       tp,
				       timeWeight,
				       DDATA(dS),
				       DDATA(x_track_ext),
				       DDATA(t_track_ext),
				       IDATA(element_track_ext),
				       IDATA(flag_track_ext),
				       IDATA(u_l2g), 
				       DDATA(u_dof),
				       DDATA(elementResidual_u), 
				       IDATA(sdInfo_u_rowptr),
				       IDATA(sdInfo_u_colind),
				       offset_u,
				       stride_u, 
				       IDATA(isFluxBoundary_u),
				       DDATA(ebqe_bc_flux_u_ext),
				       DDATA(globalResidual));

  Py_INCREF(Py_None); 
  return Py_None;
 
}




static PyMethodDef cLADR2DellamMethods[] = {
 { "calculateResidual",
    cLADR2Dellam_calculateResidual,
   METH_VARARGS, 
   "Calculate the global residual for the non-conservative level set equation"},
 { "calculateJacobian",
    cLADR2Dellam_calculateJacobian,
   METH_VARARGS, 
   "Calculate the global Jacobian for the non-conservative level set equation"},
 { "markInflowBoundaryPoints",
    cLADR2Dellam_markInflowBoundaryPoints,
   METH_VARARGS, 
   "mark points that need to be integrated for inflow boundary"},
 { "accumulateInflowFluxInGlobalResidual",
   cLADR2Dellam_accumulateInflowFluxInGlobalResidual,
   METH_VARARGS, 
   "apply inflow boundary contribution to residual "},
 { NULL,NULL,0,NULL}
};

extern "C"
{
PyMODINIT_FUNC initcLADR2Dellam(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cLADR2Dellam", cLADR2DellamMethods);
  d = PyModule_GetDict(m);
  import_array();
}
}//extern "C"
