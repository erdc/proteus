#include "Python.h"
#include "numpy/arrayobject.h"
#include "postprocessing.h"
/** \file cpostprocessingModule.h
    \defgroup cpostprocessing cpostprocessing
    \brief Python interface to velocity postprocessing library
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define NSF(p) ((NodeStarFactor *) p)

static PyObject*
cpostprocessingPostProcessRT0velocityFromP1nc(PyObject* self,
						PyObject* args)
{
  PyObject *u,*gradu,*elementBarycenters,*detJ,
    *a,*f,*r,*rt0vdofs,
    *nFreeDOF_element,*freeLocal_element,*sqrt_det_g,
    *n,*quad_a,*quad_f,*w_dV_r,
    *mt,*w_dV_m;
  mt = NULL; w_dV_m = NULL;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO|OO",
                       &nFreeDOF_element,
		       &freeLocal_element,
                       &detJ,
		       &sqrt_det_g,
		       &n,
                       &elementBarycenters,
		       &quad_a,
		       &quad_f,
		       &w_dV_r,
		       &u,
                       &gradu,
                       &a,
                       &f,
                       &r,
                       &rt0vdofs,
		       &w_dV_m,&mt))    
    return NULL;

  if (mt != NULL && w_dV_m != NULL)
    postProcessRT0velocityFromP1nc(SHAPE(gradu)[0],/*nElements_global*/
				     SHAPE(gradu)[1],/*nQuadraturePoints_element*/
				     SHAPE(w_dV_r)[2],/*nDOF_test_element*/
				     SHAPE(n)[1],    /*nElementBoundaries_element*/
				     SHAPE(n)[2],    /*nQuadraturePoints_elementBoundary*/
				     SHAPE(gradu)[2], /*spaceDim*/
				     IDATA(nFreeDOF_element),
				     IDATA(freeLocal_element),
				     DDATA(detJ),
				     DDATA(sqrt_det_g),
				     DDATA(n),
				     DDATA(elementBarycenters),
				     DDATA(quad_a),
				     DDATA(quad_f),
				     DDATA(w_dV_r),
				     DDATA(w_dV_m),
				     DDATA(u),
				     DDATA(gradu),
				     DDATA(a),
				     DDATA(f),
				     DDATA(r),
				     DDATA(mt),
				     DDATA(rt0vdofs));
  else
    postProcessRT0velocityFromP1ncNoMass(SHAPE(gradu)[0],/*nElements_global*/
					   SHAPE(gradu)[1],/*nQuadraturePoints_element*/
					   SHAPE(w_dV_r)[2],/*nDOF_test_element*/
					   SHAPE(n)[1],    /*nElementBoundaries_element*/
					   SHAPE(n)[2],    /*nQuadraturePoints_elementBoundary*/
					   SHAPE(gradu)[2], /*spaceDim*/
					   IDATA(nFreeDOF_element),
					   IDATA(freeLocal_element),
					   DDATA(detJ),
					   DDATA(sqrt_det_g),
					   DDATA(n),
					   DDATA(elementBarycenters),
					   DDATA(quad_a),
					   DDATA(quad_f),
					   DDATA(w_dV_r),
					   DDATA(u),
					   DDATA(gradu),
					   DDATA(a),
					   DDATA(f),
					   DDATA(r),
					   DDATA(rt0vdofs));
    
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingPostProcessRT0velocityFromP1nc_sd(PyObject* self,
						 PyObject* args)
{
  PyObject *u,*gradu,*elementBarycenters,*detJ,
    *a,*f,*r,*rt0vdofs,
    *nFreeDOF_element,*freeLocal_element,*sqrt_det_g,
    *n,*quad_a,*quad_f,*w_dV_r,
    *mt,*w_dV_m,*rowptr,*colind;
  mt = NULL; w_dV_m = NULL;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOO|OO",
		       &rowptr,
		       &colind,
                       &nFreeDOF_element,
		       &freeLocal_element,
                       &detJ,
		       &sqrt_det_g,
		       &n,
                       &elementBarycenters,
		       &quad_a,
		       &quad_f,
		       &w_dV_r,
		       &u,
                       &gradu,
                       &a,
                       &f,
                       &r,
                       &rt0vdofs,
		       &w_dV_m,&mt))    
    return NULL;

  if (mt != NULL && w_dV_m != NULL)
    postProcessRT0velocityFromP1nc_sd(SHAPE(gradu)[0],/*nElements_global*/
				     SHAPE(gradu)[1],/*nQuadraturePoints_element*/
				     SHAPE(w_dV_r)[2],/*nDOF_test_element*/
				     SHAPE(n)[1],    /*nElementBoundaries_element*/
				     SHAPE(n)[2],    /*nQuadraturePoints_elementBoundary*/
				     SHAPE(gradu)[2], /*spaceDim*/
				      IDATA(rowptr),
				      IDATA(colind),
				      IDATA(nFreeDOF_element),
				     IDATA(freeLocal_element),
				     DDATA(detJ),
				     DDATA(sqrt_det_g),
				     DDATA(n),
				     DDATA(elementBarycenters),
				     DDATA(quad_a),
				     DDATA(quad_f),
				     DDATA(w_dV_r),
				     DDATA(w_dV_m),
				     DDATA(u),
				     DDATA(gradu),
				     DDATA(a),
				     DDATA(f),
				     DDATA(r),
				     DDATA(mt),
				     DDATA(rt0vdofs));
  else
    postProcessRT0velocityFromP1ncNoMass_sd(SHAPE(gradu)[0],/*nElements_global*/
					   SHAPE(gradu)[1],/*nQuadraturePoints_element*/
					   SHAPE(w_dV_r)[2],/*nDOF_test_element*/
					   SHAPE(n)[1],    /*nElementBoundaries_element*/
					   SHAPE(n)[2],    /*nQuadraturePoints_elementBoundary*/
					   SHAPE(gradu)[2], /*spaceDim*/
					    IDATA(rowptr),
					    IDATA(colind),
					    IDATA(nFreeDOF_element),
					   IDATA(freeLocal_element),
					   DDATA(detJ),
					   DDATA(sqrt_det_g),
					   DDATA(n),
					   DDATA(elementBarycenters),
					   DDATA(quad_a),
					   DDATA(quad_f),
					   DDATA(w_dV_r),
					   DDATA(u),
					   DDATA(gradu),
					   DDATA(a),
					   DDATA(f),
					   DDATA(r),
					   DDATA(rt0vdofs));
    
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingUpdateRT0velocityWithAveragedPotentialP1nc(PyObject* self,
							  PyObject* args)
{
  
  PyObject *detJ, *quad_a, *phi, *gradphi, *a, *rt0vdofs_element;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &detJ,
		       &quad_a,
		       &phi,
		       &gradphi,
		       &a,
                       &rt0vdofs_element))
    return NULL;

  updateRT0velocityWithAveragedPotentialP1nc(SHAPE(a)[0],
					     SHAPE(a)[1],
					     SHAPE(gradphi)[ND(gradphi)-1],
					     DDATA(detJ),
					     DDATA(quad_a),
					     DDATA(phi),
					     DDATA(gradphi),
					     DDATA(a), 
					     DDATA(rt0vdofs_element));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingUpdateRT0velocityWithAveragedPotentialP1nc_sd(PyObject* self,
							  PyObject* args)
{
  
  PyObject *detJ, *quad_a, *phi, *gradphi, *a, *rt0vdofs_element,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
		       &rowptr,
		       &colind,
                       &detJ,
		       &quad_a,
		       &phi,
		       &gradphi,
		       &a,
                       &rt0vdofs_element))
    return NULL;

  updateRT0velocityWithAveragedPotentialP1nc_sd(SHAPE(a)[0],
						SHAPE(a)[1],
						SHAPE(gradphi)[ND(gradphi)-1],
						IDATA(rowptr),
						IDATA(colind),
						DDATA(detJ),
						DDATA(quad_a),
						DDATA(phi),
						DDATA(gradphi),
						DDATA(a), 
						DDATA(rt0vdofs_element));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetElementRT0velocityValues(PyObject* self,
					 PyObject* args)
{
  
  PyObject *x_element,*rt0vdofs_element,*v_element;
  if(!PyArg_ParseTuple(args,"OOO",
                       &x_element,
                       &rt0vdofs_element,
		       &v_element))    
    return NULL;

  getElementRT0velocityValues(SHAPE(v_element)[0],
			      SHAPE(v_element)[1],
			      SHAPE(v_element)[2],
			      DDATA(x_element),
			      DDATA(rt0vdofs_element),
			      DDATA(v_element));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cpostprocessingGetElementBoundaryRT0velocityValues(PyObject* self,
						 PyObject* args)
{
  
  PyObject *x_elementBoundary,*rt0vdofs_element,*v_elementBoundary;
  if(!PyArg_ParseTuple(args,"OOO",
                       &x_elementBoundary,
                       &rt0vdofs_element,
		       &v_elementBoundary))    
    return NULL;

  getElementBoundaryRT0velocityValues(SHAPE(v_elementBoundary)[0],
				      SHAPE(v_elementBoundary)[1],
				      SHAPE(v_elementBoundary)[2],
				      SHAPE(v_elementBoundary)[3],
				      DDATA(x_elementBoundary),
				      DDATA(rt0vdofs_element),
				      DDATA(v_elementBoundary));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetGlobalElementBoundaryRT0velocityValues(PyObject* self,
						       PyObject* args)
{
  
  PyObject *elementBoundaryElementsArray,*x_elementBoundary_global,
    *rt0vdofs_element,*v_elementBoundary_global;
  if(!PyArg_ParseTuple(args,"OOOO",
		       &elementBoundaryElementsArray,
                       &x_elementBoundary_global,
                       &rt0vdofs_element,
		       &v_elementBoundary_global))    
    return NULL;

  getGlobalElementBoundaryRT0velocityValues(SHAPE(v_elementBoundary_global)[0],
					    SHAPE(v_elementBoundary_global)[1],
					    SHAPE(v_elementBoundary_global)[2],
					    IDATA(elementBoundaryElementsArray),
					    DDATA(x_elementBoundary_global),
					    DDATA(rt0vdofs_element),
					    DDATA(v_elementBoundary_global));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cpostprocessingGetGlobalExteriorElementBoundaryRT0velocityValues(PyObject* self,
								 PyObject* args)
{
  
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,*x_elementBoundary_global,
    *rt0vdofs_element,*v_elementBoundary_global;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &exteriorElementBoundariesArray,
		       &elementBoundaryElementsArray,
                       &x_elementBoundary_global,
                       &rt0vdofs_element,
		       &v_elementBoundary_global))    
    return NULL;

  getGlobalExteriorElementBoundaryRT0velocityValues(SHAPE(exteriorElementBoundariesArray)[0],
						    SHAPE(v_elementBoundary_global)[1],
						    SHAPE(v_elementBoundary_global)[2],
						    IDATA(elementBoundaryElementsArray),
						    IDATA(exteriorElementBoundariesArray),
						    DDATA(x_elementBoundary_global),
						    DDATA(rt0vdofs_element),
						    DDATA(v_elementBoundary_global));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingPostProcessRT0potentialFromP1nc(PyObject* self,
					     PyObject* args)
{
  
  PyObject *uQuadratureWeights_element,*uQuadratureWeights_elementBoundary,
    *u,*gradu,*elementBarycenters,*aElementQuadWeights,*detJ,
    *x,*x_elementBoundary,*u_elementBoundary,
    *n,*a,*f,*r,*rt0vdofs,*rt0potential;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOO",
		       &uQuadratureWeights_element,
		       &elementBarycenters,
		       &aElementQuadWeights,
		       &detJ,
		       &uQuadratureWeights_elementBoundary,
		       &x,
		       &u,
                       &gradu,
		       &x_elementBoundary,
		       &u_elementBoundary,
		       &n,
                       &a,
                       &f,
                       &r,
                       &rt0vdofs,
		       &rt0potential))    
    return NULL;

  postProcessRT0potentialFromP1nc(SHAPE(u)[0],
				  SHAPE(u)[1],
				  SHAPE(n)[1],
				  SHAPE(n)[2],
				  SHAPE(n)[3],
				  DDATA(uQuadratureWeights_element),
				  DDATA(elementBarycenters),
				  DDATA(aElementQuadWeights),
				  DDATA(detJ),
				  DDATA(uQuadratureWeights_elementBoundary),
				  DDATA(x),
				  DDATA(u),
				  DDATA(gradu),
				  DDATA(x_elementBoundary),
				  DDATA(u_elementBoundary),
				  DDATA(n),
				  DDATA(a),
				  DDATA(f),
				  DDATA(r),
				  DDATA(rt0vdofs),
				  DDATA(rt0potential));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingPostProcessRT0potentialFromP1nc_sd(PyObject* self,
					     PyObject* args)
{
  
  PyObject *uQuadratureWeights_element,*uQuadratureWeights_elementBoundary,
    *u,*gradu,*elementBarycenters,*aElementQuadWeights,*detJ,
    *x,*x_elementBoundary,*u_elementBoundary,
    *n,*a,*f,*r,*rt0vdofs,*rt0potential,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
		       &uQuadratureWeights_element,
		       &elementBarycenters,
		       &aElementQuadWeights,
		       &detJ,
		       &uQuadratureWeights_elementBoundary,
		       &x,
		       &u,
                       &gradu,
		       &x_elementBoundary,
		       &u_elementBoundary,
		       &n,
                       &a,
                       &f,
                       &r,
                       &rt0vdofs,
		       &rt0potential))    
    return NULL;

  postProcessRT0potentialFromP1nc_sd(SHAPE(u)[0],
				     SHAPE(u)[1],
				     SHAPE(n)[1],
				     SHAPE(n)[2],
				     SHAPE(n)[3],
				     IDATA(rowptr),
				     IDATA(colind),
				     DDATA(uQuadratureWeights_element),
				     DDATA(elementBarycenters),
				     DDATA(aElementQuadWeights),
				     DDATA(detJ),
				     DDATA(uQuadratureWeights_elementBoundary),
				     DDATA(x),
				     DDATA(u),
				     DDATA(gradu),
				     DDATA(x_elementBoundary),
				     DDATA(u_elementBoundary),
				     DDATA(n),
				     DDATA(a),
				     DDATA(f),
				     DDATA(r),
				     DDATA(rt0vdofs),
				     DDATA(rt0potential));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingBuildLocalBDM1projectionMatrices(PyObject* self,
					      PyObject* args)
{
  
  PyObject *w_dS_f,*ebq_n,*ebq_v,*BDMmat_element;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &w_dS_f,
                       &ebq_n,
                       &ebq_v,
		       &BDMmat_element))
    return NULL;

  buildLocalBDM1projectionMatrices(SHAPE(ebq_n)[0],
				   SHAPE(ebq_n)[1],
				   SHAPE(ebq_n)[2],
				   SHAPE(ebq_n)[3],
				   SHAPE(w_dS_f)[3],
				   SHAPE(ebq_v)[3],
				   SHAPE(BDMmat_element)[1],
				   DDATA(w_dS_f),
				   DDATA(ebq_n),
				   DDATA(ebq_v),
				   DDATA(BDMmat_element));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingBuildLocalBDM2projectionMatrices(PyObject* self,
					      PyObject* args)
{
  
  PyObject *w_dS_f,*ebq_n,*ebq_v,*BDMmat_element,*q_basis_vals,*w_int_test_grads,*w_int_div_free,*piola_trial_fun,*degree;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
		       &degree,
                       &w_dS_f,
                       &ebq_n,
                       &ebq_v,
		       &q_basis_vals,
		       &w_int_test_grads,
		       &w_int_div_free,
		       &piola_trial_fun,
		       &BDMmat_element))
    return NULL;

  buildLocalBDM2projectionMatrices(DDATA(degree),
				   SHAPE(ebq_n)[0],
				   SHAPE(ebq_n)[1],
				   SHAPE(ebq_n)[2],
				   SHAPE(q_basis_vals)[1],
				   SHAPE(ebq_n)[3],
				   SHAPE(w_dS_f)[3],
				   SHAPE(ebq_v)[3],
				   SHAPE(w_int_test_grads)[2],
				   SHAPE(BDMmat_element)[1],
				   DDATA(w_dS_f),
				   DDATA(ebq_n),
				   DDATA(ebq_v),
				   DDATA(BDMmat_element),
				   DDATA(q_basis_vals),
				   DDATA(w_int_test_grads),
            			   DDATA(w_int_div_free),
			           DDATA(piola_trial_fun));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingFactorLocalBDM1projectionMatrices(PyObject* self,
					       PyObject* args)
{
  
  PyObject *BDMmat_element,*BDMmatPivots_element;
  if(!PyArg_ParseTuple(args,"OO",
		       &BDMmat_element,
		       &BDMmatPivots_element))
    return NULL;

  factorLocalBDM1projectionMatrices(SHAPE(BDMmat_element)[0],
				    SHAPE(BDMmat_element)[1],
				    DDATA(BDMmat_element),
				    IDATA(BDMmatPivots_element));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingFactorLocalBDM2projectionMatrices(PyObject* self,
					       PyObject* args)
{
  
  PyObject *BDMmat_element,*BDMmatPivots_element;
  if(!PyArg_ParseTuple(args,"OO",
		       &BDMmat_element,
		       &BDMmatPivots_element))
    return NULL;

  factorLocalBDM2projectionMatrices(SHAPE(BDMmat_element)[0],
				    SHAPE(BDMmat_element)[1],
				    DDATA(BDMmat_element),
				    IDATA(BDMmatPivots_element));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingSolveLocalBDM1projection(PyObject* self,
				      PyObject* args)
{
  
  PyObject *w_dS_f,*ebq_n,*ebq_velocity,*q_vdofs,*BDMmat_element,*BDMmatPivots_element;
  if(!PyArg_ParseTuple(args,"OOOOOO",
		       &BDMmat_element,
		       &BDMmatPivots_element,
                       &w_dS_f,
                       &ebq_n,
                       &ebq_velocity,
		       &q_vdofs))
    return NULL;

  solveLocalBDM1projection(SHAPE(ebq_n)[0],
			   SHAPE(ebq_n)[1],
			   SHAPE(ebq_n)[2],
			   SHAPE(ebq_n)[3],
			   SHAPE(w_dS_f)[3],
			   SHAPE(BDMmat_element)[1],
			   DDATA(BDMmat_element),
			   IDATA(BDMmatPivots_element),
			   DDATA(w_dS_f),
			   DDATA(ebq_n),
			   DDATA(ebq_velocity),
			   DDATA(q_vdofs));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingSolveLocalBDM2projection(PyObject* self,
				      PyObject* args)
{
  
  PyObject *w_dS_f,*ebq_n,*ebq_velocity,*q_vdofs,*BDMmat_element,*BDMmatPivots_element,*q_velocity,*w_interior_gradients;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
		       &BDMmat_element,
		       &BDMmatPivots_element,
                       &w_dS_f,
                       &ebq_n,
		       &w_interior_gradients,
                       &ebq_velocity,
		       &q_velocity,
		       &q_vdofs))
    return NULL;

  solveLocalBDM2projection(SHAPE(ebq_n)[0],
			   SHAPE(ebq_n)[1],
			   SHAPE(ebq_n)[2],
			   SHAPE(ebq_n)[3],
			   SHAPE(w_dS_f)[3],
			   SHAPE(BDMmat_element)[1],
			   DDATA(BDMmat_element),
			   IDATA(BDMmatPivots_element),
			   DDATA(w_dS_f),
			   DDATA(ebq_n),
			   DDATA(w_interior_gradients),
			   DDATA(q_velocity),
			   DDATA(ebq_velocity),
			   DDATA(q_vdofs));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingBuildBDM2rhs(PyObject* self,
				      PyObject* args)
{
  
  PyObject *w_dS_f,*ebq_n,*ebq_velocity,*q_vdofs,*BDMmat_element,*BDMmatPivots_element,*q_velocity,
    *w_interior_gradients, *w_interior_divfree;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
		       &BDMmat_element,
		       &BDMmatPivots_element,
                       &w_dS_f,
                       &ebq_n,
		       &w_interior_gradients,
		       &w_interior_divfree,
                       &ebq_velocity,
		       &q_velocity,
		       &q_vdofs))
    return NULL;

  buildBDM2rhs(SHAPE(ebq_n)[0],
	       SHAPE(ebq_n)[1],
	       SHAPE(ebq_n)[2],
	       SHAPE(q_velocity)[1],
	       SHAPE(ebq_n)[3],
	       SHAPE(w_dS_f)[3],
	       SHAPE(BDMmat_element)[1],
	       DDATA(BDMmat_element),
	       IDATA(BDMmatPivots_element),
	       DDATA(w_dS_f),
	       DDATA(ebq_n),
	       DDATA(w_interior_gradients),
	       DDATA(w_interior_divfree),
	       DDATA(ebq_velocity),
	       DDATA(q_velocity),		
	       DDATA(q_vdofs));

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
cpostprocessingSolveLocalBDM1projectionFromFlux(PyObject* self,
                                              PyObject* args)
{
  
  PyObject *w_dS_f,*ebq_global_flux,*q_vdofs,*BDMmat_element,*BDMmatPivots_element,*elementBoundaryElementsArray,*elementBoundariesArray;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &BDMmat_element,
		       &BDMmatPivots_element,
                       &elementBoundaryElementsArray,
                       &elementBoundariesArray,
                       &w_dS_f,
                       &ebq_global_flux,
		       &q_vdofs))
    return NULL;

  solveLocalBDM1projectionFromFlux(SHAPE(w_dS_f)[0],
                                   SHAPE(w_dS_f)[1],
                                   SHAPE(w_dS_f)[2],
                                   SHAPE(w_dS_f)[3],
                                   SHAPE(BDMmat_element)[1],
                                   DDATA(BDMmat_element),
                                   IDATA(BDMmatPivots_element),
                                   IDATA(elementBoundaryElementsArray),
                                   IDATA(elementBoundariesArray),
                                   DDATA(w_dS_f),
                                   DDATA(ebq_global_flux),
                                   DDATA(q_vdofs));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetElementBDM1velocityValuesLagrangeRep(PyObject* self,
						     PyObject* args)
{
  
  PyObject *q_v,*p1_vdofs,*q_velocity;
  int nVDOFs_element=1;
  if(!PyArg_ParseTuple(args,"OOO",
		       &q_v,
		       &p1_vdofs,
		       &q_velocity))
    return NULL;
  nVDOFs_element = SHAPE(q_velocity)[2]*(SHAPE(q_velocity)[2]+1);
  
  if (ND(p1_vdofs) > 1)
    assert(nVDOFs_element == SHAPE(p1_vdofs)[1]);

  getElementBDM1velocityValuesLagrangeRep(SHAPE(q_velocity)[0],
					  SHAPE(q_velocity)[1],
					  SHAPE(q_velocity)[2],
					  SHAPE(q_v)[2],
					  nVDOFs_element,
					  DDATA(q_v),
					  DDATA(p1_vdofs),
					  DDATA(q_velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetElementBDM2velocityValuesLagrangeRep(PyObject* self,
						     PyObject* args)
{
  
  PyObject *q_v,*p1_vdofs,*q_velocity;
  int nVDOFs_element=1;
  if(!PyArg_ParseTuple(args,"OOO",
		       &q_v,
		       &p1_vdofs,
		       &q_velocity))
    return NULL;
  // this is really the dimension ... there should be a better way
  // to handle this ...
  nVDOFs_element = (SHAPE(q_velocity)[2]+1)*(SHAPE(q_velocity)[2]+2);
  
  if (ND(p1_vdofs) > 1)
    assert(nVDOFs_element == SHAPE(p1_vdofs)[1]);

  getElementBDM2velocityValuesLagrangeRep(SHAPE(q_velocity)[0],
					  SHAPE(q_velocity)[1],
					  SHAPE(q_velocity)[2],
					  SHAPE(q_v)[2],
					  nVDOFs_element,
					  DDATA(q_v),
					  DDATA(p1_vdofs),
					  DDATA(q_velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

// tjp added for LDG velocity local representation for error analysis
static PyObject*
cpostprocessingGetElementLDGvelocityValuesLagrangeRep(PyObject* self,
						     PyObject* args)
{
  
  PyObject *q_v,*vdofs,*q_velocity;
  int nVDOF_element=1;
  int nDOF_trial_element=1;
  
  if(!PyArg_ParseTuple(args,"OOO",
		       &q_v,
		       &vdofs,
		       &q_velocity))
    return NULL;
  
  nVDOF_element =SHAPE(vdofs)[1];
  nDOF_trial_element=SHAPE(q_v)[2];
  //printf('nVDOFs_element=%d',nVDOFs_element);
  //nVDOFs_element = SHAPE(q_velocity)[2]*(SHAPE(q_velocity)[2]+1);
  if (ND(vdofs) > 1)
    assert(nVDOF_element == SHAPE(q_v)[2]*SHAPE(q_velocity)[2]);

  getElementLDGvelocityValuesLagrangeRep(SHAPE(q_velocity)[0],
					  SHAPE(q_velocity)[1],
					  SHAPE(q_velocity)[2],
					  nDOF_trial_element,
					  nVDOF_element,
					  DDATA(q_v),
					  DDATA(vdofs),
					  DDATA(q_velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(PyObject* self,
									     PyObject* args)
{
  
  PyObject *elementBoundaryElementsArray,
    *exteriorElementBoundariesArray,
    *ebqe_v,*p1_vdofs,*ebqe_velocity;
  int nVDOFs_element=1;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &elementBoundaryElementsArray,
		       &exteriorElementBoundariesArray,
		       &ebqe_v,
		       &p1_vdofs,
		       &ebqe_velocity))
    return NULL;
  nVDOFs_element = SHAPE(ebqe_velocity)[2]*(SHAPE(ebqe_velocity)[2]+1);
  if (ND(p1_vdofs) > 1)
    assert(nVDOFs_element == SHAPE(p1_vdofs)[1]);

  getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(SHAPE(ebqe_velocity)[0],
								SHAPE(ebqe_velocity)[1],
								SHAPE(ebqe_velocity)[2],
								SHAPE(ebqe_v)[2],
								nVDOFs_element,
								IDATA(elementBoundaryElementsArray),
								IDATA(exteriorElementBoundariesArray),
								DDATA(ebqe_v),
								DDATA(p1_vdofs),
								DDATA(ebqe_velocity));
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
cpostprocessingGetGlobalElementBoundaryBDM1velocityValuesLagrangeRep(PyObject* self,
									     PyObject* args)
{
  
  PyObject *elementBoundaryElementsArray,
    *exteriorElementBoundariesArray,
    *ebqe_v,*p1_vdofs,*ebq_global_velocity;
  int nVDOFs_element=1;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &elementBoundaryElementsArray,
		       &exteriorElementBoundariesArray,
		       &ebqe_v,
		       &p1_vdofs,
		       &ebq_global_velocity))
    return NULL;
  nVDOFs_element = SHAPE(ebq_global_velocity)[2]*(SHAPE(ebq_global_velocity)[2]+1);
  if (ND(p1_vdofs) > 1)
    assert(nVDOFs_element == SHAPE(p1_vdofs)[1]);

  getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep(SHAPE(ebqe_v)[0],
								SHAPE(ebq_global_velocity)[1],
								SHAPE(ebq_global_velocity)[2],
								SHAPE(ebqe_v)[2],
								nVDOFs_element,
								IDATA(elementBoundaryElementsArray),
								IDATA(exteriorElementBoundariesArray),
								DDATA(ebqe_v),
								DDATA(p1_vdofs),
								DDATA(ebq_global_velocity));
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject*
cpostprocessingGetElementBoundaryBDM1velocityValuesLagrangeRep(PyObject* self,
									     PyObject* args)
{
  
  PyObject *elementBoundaryElementsArray,
    *exteriorElementBoundariesArray,
    *ebq_v,*p1_vdofs,*ebq_velocity;
  int nVDOFs_element=1;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &elementBoundaryElementsArray,
		       &exteriorElementBoundariesArray,
		       &ebq_v,
		       &p1_vdofs,
		       &ebq_velocity))
    return NULL;
  nVDOFs_element = SHAPE(ebq_velocity)[3]*(SHAPE(ebq_velocity)[3]+1);
  if (ND(p1_vdofs) > 1)
    assert(nVDOFs_element == SHAPE(p1_vdofs)[1]);

  getElementBoundaryBDM1velocityValuesLagrangeRep(SHAPE(ebq_velocity)[0],
								SHAPE(ebq_velocity)[1],
								SHAPE(ebq_velocity)[2],
								SHAPE(ebq_velocity)[3],
								SHAPE(ebq_v)[2],
								nVDOFs_element,
								IDATA(elementBoundaryElementsArray),
								IDATA(exteriorElementBoundariesArray), /*need the correct mapping here */
								DDATA(ebq_v),
								DDATA(p1_vdofs),
								DDATA(ebq_velocity));
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
cpostprocessingProjectElementBoundaryVelocityToRT0fluxRep(PyObject* self,
							PyObject* args)
{

  
  PyObject *elementBoundaryQuadratureWeights,*n,*v_elementBoundary,
    *rt0vdofs_element;
  if(!PyArg_ParseTuple(args,"OOOO",
		       &elementBoundaryQuadratureWeights,
		       &n,
		       &v_elementBoundary,
		       &rt0vdofs_element))
    return NULL;

  projectElementBoundaryVelocityToRT0fluxRep(SHAPE(n)[0],
					     SHAPE(n)[1],
					     SHAPE(n)[2],
					     SHAPE(n)[3],
					     DDATA(elementBoundaryQuadratureWeights),
					     DDATA(n),
					     DDATA(v_elementBoundary),
					     DDATA(rt0vdofs_element));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingProjectElementBoundaryFluxToRT0fluxRep(PyObject* self,
                                                    PyObject* args)
{
  
  PyObject *elementBoundaryElementsArray,
    *elementBoundariesArray,
    *elementBoundaryQuadratureWeights,
    *flux_elementBoundary,
    *rt0vdofs_element;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &elementBoundaryElementsArray,
                       &elementBoundariesArray,
		       &elementBoundaryQuadratureWeights,
		       &flux_elementBoundary,
		       &rt0vdofs_element))
    return NULL;
  projectElementBoundaryFluxToRT0fluxRep(SHAPE(elementBoundaryQuadratureWeights)[0],
                                         SHAPE(elementBoundaryQuadratureWeights)[1],
                                         SHAPE(elementBoundaryQuadratureWeights)[2],
                                         SHAPE(rt0vdofs_element)[1],
                                         IDATA(elementBoundaryElementsArray),
                                         IDATA(elementBoundariesArray),
                                         DDATA(elementBoundaryQuadratureWeights),
                                         DDATA(flux_elementBoundary),
                                         DDATA(rt0vdofs_element));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetElementRT0velocityValuesFluxRep(PyObject* self,
						PyObject* args)
{
  
  PyObject *nodeArray,*elementNodesArray,*abs_det_J,*x_element,
    *rt0vdofs_element,*v_element;
  if(!PyArg_ParseTuple(args,"OOOOOO",
		       &nodeArray,
		       &elementNodesArray,
		       &abs_det_J,
		       &x_element,
		       &rt0vdofs_element,
		       &v_element))
    return NULL;

  getElementRT0velocityValuesFluxRep(SHAPE(v_element)[0],
				     SHAPE(elementNodesArray)[1],
				     SHAPE(v_element)[1],
				     SHAPE(v_element)[2],
				     SHAPE(abs_det_J)[1],
				     DDATA(nodeArray),
				     IDATA(elementNodesArray),
				     DDATA(abs_det_J),
				     DDATA(x_element),
				     DDATA(rt0vdofs_element),
				     DDATA(v_element));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetElementBoundaryRT0velocityValuesFluxRep(PyObject* self,
                                                        PyObject* args)
{
  
  PyObject *nodeArray,
    *elementNodesArray,
    *abs_det_J,
    *x_elementBoundary,
    *rt0vdofs_element,
    *v_elementBoundary;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &nodeArray,
                       &elementNodesArray,
                       &abs_det_J,
                       &x_elementBoundary,
                       &rt0vdofs_element,
		       &v_elementBoundary))    
    return NULL;
  getElementBoundaryRT0velocityValuesFluxRep(SHAPE(v_elementBoundary)[0],
                                             SHAPE(v_elementBoundary)[1],
                                             SHAPE(v_elementBoundary)[2],
                                             SHAPE(v_elementBoundary)[3],
                                             SHAPE(abs_det_J)[1],
                                             DDATA(nodeArray),
                                             IDATA(elementNodesArray),
                                             DDATA(abs_det_J),
                                             DDATA(x_elementBoundary),
                                             DDATA(rt0vdofs_element),
                                             DDATA(v_elementBoundary));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetGlobalElementBoundaryRT0velocityValuesFluxRep(PyObject* self,
                                                              PyObject* args)
{
  
  PyObject *nodeArray,
    *elementNodesArray,
    *abs_det_J,
    *elementBoundaryElementsArray,
    *x_elementBoundary_global,
    *rt0vdofs_element,
    *v_elementBoundary_global;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &nodeArray,
                       &elementNodesArray,
                       &elementBoundaryElementsArray,
                       &abs_det_J,
                       &x_elementBoundary_global,
                       &rt0vdofs_element,
		       &v_elementBoundary_global))    
    return NULL;

  getGlobalElementBoundaryRT0velocityValuesFluxRep(SHAPE(v_elementBoundary_global)[0],
                                                   SHAPE(v_elementBoundary_global)[1],
                                                   SHAPE(v_elementBoundary_global)[2],
                                                   SHAPE(abs_det_J)[1],
                                                   DDATA(nodeArray),
                                                   IDATA(elementNodesArray),
                                                   IDATA(elementBoundaryElementsArray),
                                                   DDATA(abs_det_J),
                                                   DDATA(x_elementBoundary_global),
                                                   DDATA(rt0vdofs_element),
                                                   DDATA(v_elementBoundary_global));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cpostprocessingGetGlobalExteriorElementBoundaryRT0velocityValuesFluxRep(PyObject* self,
									PyObject* args)
{
  
  PyObject *nodeArray,
    *elementNodesArray,
    *abs_det_J,
    *elementBoundaryElementsArray,
    *exteriorElementBoundariesArray,
    *x_ebqe,
    *rt0vdofs_element,
    *v_ebqe;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
		       &nodeArray,
                       &elementNodesArray,
		       &elementBoundaryElementsArray,
                       &exteriorElementBoundariesArray,
                       &abs_det_J,
                       &x_ebqe,
                       &rt0vdofs_element,
		       &v_ebqe))
    return NULL;

  getGlobalExteriorElementBoundaryRT0velocityValuesFluxRep(SHAPE(v_ebqe)[0],
							   SHAPE(v_ebqe)[1],
							   SHAPE(v_ebqe)[2],
							   SHAPE(abs_det_J)[1],
							   DDATA(nodeArray),
							   IDATA(elementNodesArray),
							   IDATA(elementBoundaryElementsArray),
							   IDATA(exteriorElementBoundariesArray),
							   DDATA(abs_det_J),
							   DDATA(x_ebqe),
							   DDATA(rt0vdofs_element),
							   DDATA(v_ebqe));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingGetRT0velocityValuesFluxRep_arbitraryElementMembership(PyObject* self,
								      PyObject* args)
{
  
  PyObject *nodeArray,*elementNodesArray,*abs_det_J,*x,
    *element_locations,*rt0vdofs_element,*v_element;
  int i,nPoints=1;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &nodeArray,
		       &elementNodesArray,
		       &abs_det_J,
		       &x,
		       &element_locations,
		       &rt0vdofs_element,
		       &v_element))
    return NULL;
  for (i=0; i < ND(x)-1; i++)
    nPoints *= SHAPE(x)[i];
  getRT0velocityValuesFluxRep_arbitraryElementMembership(SHAPE(elementNodesArray)[0],
							 SHAPE(elementNodesArray)[1],
							 nPoints,
							 SHAPE(v_element)[1],
							 SHAPE(abs_det_J)[1],
							 DDATA(nodeArray),
							 IDATA(elementNodesArray),
							 DDATA(abs_det_J),
							 DDATA(x),
							 IDATA(element_locations),
							 DDATA(rt0vdofs_element),
							 DDATA(v_element));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingFluxCorrectionVelocityUpdate(PyObject* self, 
				      PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* dS;
  PyObject* n;
  PyObject* fluxCorrection;
  PyObject* velocity;
  PyObject* velocity_element;

  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &dS,
		       &n,
		       &fluxCorrection,
		       &velocity,
		       &velocity_element))

    return NULL;
  fluxCorrectionVelocityUpdate(SHAPE(dS)[0],
			       SHAPE(velocity)[0],
			       SHAPE(interiorElementBoundaries)[0],
			       SHAPE(exteriorElementBoundaries)[0],
			       SHAPE(dS)[1],
			       SHAPE(dS)[2],
			       SHAPE(velocity)[2],
			       IDATA(interiorElementBoundaries),
			       IDATA(exteriorElementBoundaries),
			       IDATA(elementBoundaryElements),
			       IDATA(elementBoundaryLocalElementBoundaries),
			       DDATA(dS),
			       DDATA(n),
			       DDATA(fluxCorrection),
			       DDATA(velocity),
			       DDATA(velocity_element));


  Py_INCREF(Py_None); 
  return Py_None;
}



static PyObject*
cpostprocessingComputeFluxCorrectionPWC(PyObject* self, 
				  PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* pwcV;
  PyObject* pwcW;
  PyObject* fluxCorrection;

  if(!PyArg_ParseTuple(args,"OOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &pwcW,
		       &pwcV,
		       &fluxCorrection))

    return NULL;
  computeFluxCorrectionPWC(SHAPE(elementBoundaryElements)[0],
			   SHAPE(interiorElementBoundaries)[0],
			   SHAPE(exteriorElementBoundaries)[0],
			   IDATA(interiorElementBoundaries),
			   IDATA(exteriorElementBoundaries),
			   IDATA(elementBoundaryElements),
			   DDATA(pwcW),
			   DDATA(pwcV),
			   DDATA(fluxCorrection));


  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingSunWheelerGSsweep(PyObject* self, 
			   PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* dS;
  PyObject* n;
  PyObject* sqrt_det_g;
  PyObject* alphaFactor;
  PyObject* fluxCorrection;
  PyObject* conservationResidual;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &dS,
		       &n,
		       &sqrt_det_g,
		       &alphaFactor,
		       &fluxCorrection,
		       &conservationResidual))
    return NULL;
  sunWheelerGSsweep(SHAPE(conservationResidual)[0],
		    SHAPE(elementBoundaryElements)[0],
		    SHAPE(interiorElementBoundaries)[0],
		    SHAPE(exteriorElementBoundaries)[0],
		    SHAPE(dS)[1],
		    SHAPE(dS)[2],
		    SHAPE(n)[2],
		    IDATA(interiorElementBoundaries),
		    IDATA(exteriorElementBoundaries),
		    IDATA(elementBoundaryElements),
		    IDATA(elementBoundaryLocalElementBoundaries),
		    DDATA(dS),
		    DDATA(n),
		    DDATA(sqrt_det_g),
		    DDATA(alphaFactor),
		    DDATA(fluxCorrection),
		    DDATA(conservationResidual));

  Py_INCREF(Py_None); 
  return Py_None;
}

/***********************************************************************
  put in data structure for node star solves now
 ***********************************************************************/
static PyObject*
NodeStarFactor_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  NodeStarFactor *self;
  self = (NodeStarFactor *)type->tp_alloc(type,0);
  return (PyObject*)self;
}
static int
NodeStarFactor_init(NodeStarFactor* self, PyObject *args, PyObject *kwds)
{
  PyObject *nElements_node,*nodeStarElementsArray,*nodeStarElementsNeighborsArray;
  if(!PyArg_ParseTuple(args,
                       "OOO",
                       &nElements_node,
		       &nodeStarElementsArray,
		       &nodeStarElementsNeighborsArray))
    return -1;
  
  if (nodeStar_init(SHAPE(nodeStarElementsArray)[0],
		    SHAPE(nodeStarElementsArray)[1],
		    SHAPE(nElements_node)[0],
		    IDATA(nElements_node),
		    IDATA(nodeStarElementsArray),
		    IDATA(nodeStarElementsNeighborsArray),
		    &self->N,
		    &self->subdomain_dim,
		    &self->subdomain_L,
		    &self->subdomain_R,
		    &self->subdomain_U,
		    &self->subdomain_pivots,
		    &self->subdomain_column_pivots))
    {
      PyErr_NoMemory();
      return -1;
    }
  return 0;
}
static  void
NodeStarFactor_dealloc(NodeStarFactor* self)
{
  nodeStar_free(self->N, 
		self->subdomain_dim,
		self->subdomain_L,
		self->subdomain_R,
		self->subdomain_U,
		self->subdomain_pivots,
		self->subdomain_column_pivots);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject*
NodeStarFactor_setU(NodeStarFactor* self, PyObject* args)
{

  double val;
  if (!PyArg_ParseTuple(args,
			"d",
			&val))
    return NULL;
  assert(self);
  nodeStar_setU(self,val);
  
  Py_INCREF(Py_None); 
  return Py_None;

}

static PyObject*
NodeStarFactor_copyData(NodeStarFactor* self, PyObject* args)
{

  PyObject *other;
  if (!PyArg_ParseTuple(args,
			"O",
			&other))
    return NULL;
  assert(self);
  assert(other);

  nodeStar_copy(NSF(other)->N, 
		NSF(other)->subdomain_dim,
		NSF(other)->subdomain_L,
		NSF(other)->subdomain_R,
		NSF(other)->subdomain_U,
		NSF(other)->subdomain_pivots,
		NSF(other)->subdomain_column_pivots,
		&self->N, 
		&self->subdomain_dim,
		&self->subdomain_L,
		&self->subdomain_R,
		&self->subdomain_U,
		&self->subdomain_pivots,
		&self->subdomain_column_pivots);
 		
  
  Py_INCREF(Py_None); 
  return Py_None;

}

static PyMethodDef NodeStarFactorMethods[] = {
  {"setU",
   (PyCFunction)NodeStarFactor_setU,
   METH_VARARGS,
   "set subdomain_U to constant value"},
  {"copyData",
   (PyCFunction)NodeStarFactor_copyData,
   METH_VARARGS,
   "copy other subdomain's data, only reallocate if sizes different"},
  {NULL} /*end*/
};

static PyTypeObject NodeStarFactorType = {    
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "nodeStarFactor.NodeStarFactor",             /*tp_name*/
  sizeof(NodeStarFactor), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)NodeStarFactor_dealloc,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,        /*tp_flags*/
  "node star factorization objects",           /* tp_doc */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  NodeStarFactorMethods,     /* tp_methods */
  0,                         /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)NodeStarFactor_init,      /* tp_init */
  0,                         /* tp_alloc */
  NodeStarFactor_new,                 /* tp_new */
};

/***********************************************************************
 end node star factor type
 ***********************************************************************/

static PyObject*
cpostprocessingCalculateConservationResidualPWL(PyObject* self, 
					    PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* dofMapl2g;
  PyObject* nodeStarElements;
  PyObject* nodeStarElementNeighbors;
  PyObject* nElements_node;
  PyObject* elementResidual;
  PyObject* vAverage;
  PyObject* dx;
  PyObject* w;
  PyObject* n;
  PyObject* conservationResidual;
  PyObject* vConservative;
  PyObject* vConservative_element;
  PyObject* fluxElementBoundaries;
  PyObject* nodeStarFactor;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
    		       &dofMapl2g,
		       &nodeStarElements,
		       &nodeStarElementNeighbors,
		       &nElements_node,
		       &fluxElementBoundaries,
		       &elementResidual,
		       &vAverage,
		       &dx,
		       &w,
		       &n,
		       &nodeStarFactor,
		       &conservationResidual,
		       &vConservative,
                       &vConservative_element))
    
    return NULL;
  calculateConservationResidualPWL(SHAPE(w)[0] 
				   /* nElements_global*/, 
				   SHAPE(interiorElementBoundaries)[0] 
				   /*nInteriorElementBoundaries_global*/,
				   SHAPE(exteriorElementBoundaries)[0]
				   /*nExteriorElementBoundaries_global*/,
				   SHAPE(w)[1]
				   /*nElementBoundaries_element*/,
				   SHAPE(w)[2]
				   /*nQuadraturePoints_elementBoundary*/,
				   SHAPE(w)[3]
				   /*nNodes_element*/,
				   SHAPE(dofMapl2g)[1]
				   /*nDOF_element*/,
				   SHAPE(n)[2]
				   /*nSpace */,
				   IDATA(interiorElementBoundaries),
				   IDATA(exteriorElementBoundaries),
				   IDATA(elementBoundaryElements),
				   IDATA(elementBoundaryLocalElementBoundaries),
				   IDATA(elementNodes),
				   IDATA(dofMapl2g),
				   IDATA(nodeStarElements),
				   IDATA(nodeStarElementNeighbors),
				   IDATA(nElements_node),
				   IDATA(fluxElementBoundaries),
				   DDATA(elementResidual),
				   DDATA(vAverage),
				   DDATA(dx),
				   DDATA(w),
				   DDATA(n),
				   NSF(nodeStarFactor),
				   DDATA(conservationResidual),
				   DDATA(vConservative),
				   DDATA(vConservative_element));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationJacobianPWL(PyObject* self, 
					    PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* dofMapl2g;
  PyObject* dofStarElements;
  PyObject* dofStarElementNeighbors;
  PyObject* nElements_dof;
  PyObject* internalNodes;
  PyObject* nodeStarFactor;
  PyObject* w;
  PyObject* n;
  PyObject* fluxElementBoundaries;
  PyObject* fluxBoundaryNodes;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
		       &dofMapl2g,
		       &dofStarElements,
		       &dofStarElementNeighbors,
		       &nElements_dof,
                       &internalNodes,
		       &fluxElementBoundaries,
		       &fluxBoundaryNodes,
		       &w,
		       &n,
		       &nodeStarFactor))
    
    return NULL;
  calculateConservationJacobianPWL(SHAPE(nElements_dof)[0],
				     SHAPE(internalNodes)[0],
				     SHAPE(w)[0],
				     SHAPE(interiorElementBoundaries)[0],
				     SHAPE(exteriorElementBoundaries)[0],
				     SHAPE(w)[1],
				     SHAPE(w)[2],
				     SHAPE(w)[3],
				     SHAPE(dofMapl2g)[1],
				     SHAPE(n)[2],
				     IDATA(interiorElementBoundaries),
				     IDATA(exteriorElementBoundaries),
				     IDATA(elementBoundaryElements),
				     IDATA(elementBoundaryLocalElementBoundaries),
				     IDATA(elementNodes),
  				     IDATA(dofMapl2g),
				     IDATA(dofStarElements),
				     IDATA(dofStarElementNeighbors),
				     IDATA(nElements_dof),
				     IDATA(internalNodes),
				     IDATA(fluxElementBoundaries),
				     IDATA(fluxBoundaryNodes),
				     DDATA(w),
				     DDATA(n),
				     NSF(nodeStarFactor));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationFluxPWL(PyObject* self, 
					PyObject* args)
{
  PyObject *nElements_node,
    *internalNodes,
    *nodeStarFactor;
  PyObject *fluxBoundaryNodes;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &nElements_node,
                       &internalNodes,
		       &fluxBoundaryNodes,
                       &nodeStarFactor))
    return NULL;
  calculateConservationFluxPWL(SHAPE(nElements_node)[0],
                               SHAPE(internalNodes)[0],
                               IDATA(nElements_node),
                               IDATA(internalNodes),
                               IDATA(fluxBoundaryNodes),
                               NSF(nodeStarFactor));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject*
cpostprocessingCalculateConservationFluxPWL_noNeumannFix(PyObject* self, 
							 PyObject* args)
{
  PyObject *nElements_node,
    *nodeStarFactor;
  if(!PyArg_ParseTuple(args,"OO",
                       &nElements_node,
                       &nodeStarFactor))
    return NULL;
  calculateConservationFluxPWL_noNeumannFix(SHAPE(nElements_node)[0],
					    IDATA(nElements_node),
					    NSF(nodeStarFactor));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationResidualPWL_opt(PyObject* self, 
						    PyObject* args)
{
  int nNodes_owned;
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* nodeStarElements;
  PyObject* nodeStarElementNeighbors;
  PyObject* nElements_node;
  PyObject* elementResidual;
  PyObject* vAverage;
  PyObject* dx;
  PyObject* w;
  PyObject* n;
  PyObject* conservationResidual;
  PyObject* vConservative;
  PyObject* vConservative_element;
  PyObject* fluxElementBoundaries;
  PyObject* nodeStarFactor;
  if(!PyArg_ParseTuple(args,"iOOOOOOOOOOOOOOOOOO",
		       &nNodes_owned,
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
		       &nodeStarElements,
		       &nodeStarElementNeighbors,
		       &nElements_node,
		       &fluxElementBoundaries,
		       &elementResidual,
		       &vAverage,
		       &dx,
		       &w,
		       &n,
		       &nodeStarFactor,
		       &conservationResidual,
		       &vConservative,
                       &vConservative_element))
    
    return NULL;
  calculateConservationResidualPWL_opt(nNodes_owned,
				   SHAPE(w)[0],
				   SHAPE(interiorElementBoundaries)[0],
				   SHAPE(exteriorElementBoundaries)[0],
				   SHAPE(w)[1],
				   SHAPE(w)[2],
				   SHAPE(w)[3],
				   SHAPE(n)[2],
				   IDATA(interiorElementBoundaries),
				   IDATA(exteriorElementBoundaries),
				   IDATA(elementBoundaryElements),
				   IDATA(elementBoundaryLocalElementBoundaries),
				   IDATA(elementNodes),
				   IDATA(nodeStarElements),
				   IDATA(nodeStarElementNeighbors),
				   IDATA(nElements_node),
				   IDATA(fluxElementBoundaries),
				   DDATA(elementResidual),
				   DDATA(vAverage),
				   DDATA(dx),
				   DDATA(w),
				   DDATA(n),
				   NSF(nodeStarFactor),
				   DDATA(conservationResidual),
				   DDATA(vConservative),
				   DDATA(vConservative_element));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationResidualPWL_primative(PyObject* self, 
							  PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* skipflag_elementBoundaries;
  PyObject* elementResidual;
  PyObject* dx;
  PyObject* n;
  PyObject* conservationResidual;
  PyObject* vConservative;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &skipflag_elementBoundaries,
		       &elementResidual,
		       &dx,
		       &n,
		       &conservationResidual,
		       &vConservative))
    
    return NULL;
  calculateConservationResidualPWL_primative(SHAPE(dx)[0],
					     SHAPE(interiorElementBoundaries)[0],
					     SHAPE(exteriorElementBoundaries)[0],
					     SHAPE(dx)[1],
					     SHAPE(dx)[2],
					     SHAPE(elementResidual)[1],
					     SHAPE(n)[2],
					     IDATA(interiorElementBoundaries),
					     IDATA(exteriorElementBoundaries),
					     IDATA(elementBoundaryElements),
					     IDATA(elementBoundaryLocalElementBoundaries),
					     IDATA(skipflag_elementBoundaries),
					     DDATA(elementResidual),
					     DDATA(dx),
					     DDATA(n),
					     DDATA(conservationResidual),
					     DDATA(vConservative));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationJacobianPWL_opt(PyObject* self, 
					    PyObject* args)
{
  int nNodes_owned;
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* nodeStarElements;
  PyObject* nodeStarElementNeighbors;
  PyObject* nElements_node;
  PyObject* internalNodes;
  PyObject* nodeStarFactor;
  PyObject* w;
  PyObject* n;
  PyObject* fluxElementBoundaries;
  PyObject* fluxBoundaryNodes;

  if(!PyArg_ParseTuple(args,"iOOOOOOOOOOOOOO",
		       &nNodes_owned,
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
		       &nodeStarElements,
		       &nodeStarElementNeighbors,
		       &nElements_node,
                       &internalNodes,
		       &fluxElementBoundaries,
		       &fluxBoundaryNodes,
		       &w,
		       &n,
		       &nodeStarFactor))
    
    return NULL;
  calculateConservationJacobianPWL_opt(nNodes_owned,
				   SHAPE(nElements_node)[0],
				     SHAPE(internalNodes)[0],
				     SHAPE(w)[0],
				     SHAPE(interiorElementBoundaries)[0],
				     SHAPE(exteriorElementBoundaries)[0],
				     SHAPE(w)[1],
				     SHAPE(w)[2],
				     SHAPE(w)[3],
				     SHAPE(n)[2],
				     IDATA(interiorElementBoundaries),
				     IDATA(exteriorElementBoundaries),
				     IDATA(elementBoundaryElements),
				     IDATA(elementBoundaryLocalElementBoundaries),
				     IDATA(elementNodes),
				     IDATA(nodeStarElements),
				     IDATA(nodeStarElementNeighbors),
				     IDATA(nElements_node),
				     IDATA(internalNodes),
				     IDATA(fluxElementBoundaries),
				     IDATA(fluxBoundaryNodes),
				     DDATA(w),
				     DDATA(n),
				     NSF(nodeStarFactor));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationFluxPWL_opt(PyObject* self, 
					PyObject* args)
{
  int nNodes_owned;
  PyObject *nElements_node,
    *internalNodes,
    *nodeStarFactor;
  PyObject *fluxBoundaryNodes;
  if(!PyArg_ParseTuple(args,"iOOOO",
		       &nNodes_owned,
                       &nElements_node,
                       &internalNodes,
		       &fluxBoundaryNodes,
                       &nodeStarFactor))
    return NULL;
  calculateConservationFluxPWL_opt(nNodes_owned,
			       SHAPE(nElements_node)[0],
                               SHAPE(internalNodes)[0],
                               IDATA(nElements_node),
                               IDATA(internalNodes),
                               IDATA(fluxBoundaryNodes),
                               NSF(nodeStarFactor));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationResidualPWL_interiorBoundaries(PyObject* self, 
								   PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* nodeStarElements;
  PyObject* nodeStarElementNeighbors;
  PyObject* nElements_node;
  PyObject* elementResidual;
  PyObject* vAverage;
  PyObject* dx;
  PyObject* w;
  PyObject* n;
  PyObject* conservationResidual;
  PyObject* vConservative;
  PyObject* vConservative_element;
  PyObject* fluxElementBoundaries;
  PyObject* nodeStarFactor;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
		       &nodeStarElements,
		       &nodeStarElementNeighbors,
		       &nElements_node,
		       &fluxElementBoundaries,
		       &elementResidual,
		       &vAverage,
		       &dx,
		       &w,
		       &n,
		       &nodeStarFactor,
		       &conservationResidual,
		       &vConservative,
                       &vConservative_element))
    
    return NULL;
  assert(SHAPE(fluxElementBoundaries)[0] == SHAPE(elementBoundaryLocalElementBoundaries)[0]);

  calculateConservationResidualPWL_interiorBoundaries(SHAPE(w)[0],
						      SHAPE(interiorElementBoundaries)[0],
						      SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(w)[1],
						      SHAPE(w)[2],
						      SHAPE(w)[3],
						      SHAPE(n)[2],
						      IDATA(interiorElementBoundaries),
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(elementNodes),
						      IDATA(nodeStarElements),
						      IDATA(nodeStarElementNeighbors),
						      IDATA(nElements_node),
						      IDATA(fluxElementBoundaries),
						      DDATA(elementResidual),
						      DDATA(vAverage),
						      DDATA(dx),
						      DDATA(w),
						      DDATA(n),
						      NSF(nodeStarFactor),
						      DDATA(conservationResidual),
						      DDATA(vConservative),
						      DDATA(vConservative_element));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationJacobianPWL_interiorBoundaries(PyObject* self, 
								   PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* elementNodes;
  PyObject* nodeStarElements;
  PyObject* nodeStarElementNeighbors;
  PyObject* nElements_node;
  PyObject* nodeStarFactor;
  PyObject* w;
  PyObject* n;
  PyObject* fluxElementBoundaries;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &elementNodes,
		       &nodeStarElements,
		       &nodeStarElementNeighbors,
		       &nElements_node,
		       &fluxElementBoundaries,
		       &w,
		       &n,
		       &nodeStarFactor))
    
    return NULL;

  assert(SHAPE(fluxElementBoundaries)[0] == SHAPE(elementBoundaryLocalElementBoundaries)[0]);

  calculateConservationJacobianPWL_interiorBoundaries(SHAPE(nElements_node)[0],
						      SHAPE(w)[0],
						      SHAPE(interiorElementBoundaries)[0],
						      SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(w)[1],
						      SHAPE(w)[2],
						      SHAPE(w)[3],
						      SHAPE(n)[2],
						      IDATA(interiorElementBoundaries),
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(elementNodes),
						      IDATA(nodeStarElements),
						      IDATA(nodeStarElementNeighbors),
						      IDATA(nElements_node),
						      IDATA(fluxElementBoundaries),
						      DDATA(w),
						      DDATA(n),
						      NSF(nodeStarFactor));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessing_subdomain_U_copy_global2local(PyObject* self, 
                                              PyObject* args)
{
  int max_nN_owned;
  PyObject *elementNodes,*nodeStarElements,*nodeStarFactor,*subdomain_U;
  if(!PyArg_ParseTuple(args,"iOOOO",
                       &max_nN_owned,
                       &elementNodes,
                       &nodeStarElements,
                       &nodeStarFactor,
                       &subdomain_U))
    return NULL;
  subdomain_U_copy_global2local(max_nN_owned,
                                SHAPE(subdomain_U)[0],
                                SHAPE(subdomain_U)[1],
                                IDATA(elementNodes),
                                IDATA(nodeStarElements),
                                NSF(nodeStarFactor),
                                DDATA(subdomain_U));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessing_subdomain_U_copy_local2global(PyObject* self, 
                                              PyObject* args)
{
  int max_nN_owned;
  PyObject *elementNodes,*nodeStarElements,*nodeStarFactor,*subdomain_U;
  if(!PyArg_ParseTuple(args,"iOOOO",
                       &max_nN_owned,
                       &elementNodes,
                       &nodeStarElements,
                       &nodeStarFactor,
                       &subdomain_U))
    return NULL;
  subdomain_U_copy_local2global(max_nN_owned,
                                SHAPE(subdomain_U)[0],
                                SHAPE(subdomain_U)[1],
                                IDATA(elementNodes),
                                IDATA(nodeStarElements),
                                NSF(nodeStarFactor),
                                DDATA(subdomain_U));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingCalculateConservationResidualGlobalBoundaries(PyObject* self, 
							     PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* exteriorElementBoundariesToSkip;
  PyObject* dS;
  PyObject* n;
  PyObject* elementResidual;
  PyObject* velocity;
  PyObject* conservationResidual;
 
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &exteriorElementBoundariesToSkip,
		       &dS,
		       &n,
		       &elementResidual,
		       &velocity,
		       &conservationResidual))
    return NULL;
  calculateConservationResidualGlobalBoundaries(SHAPE(conservationResidual)[0],
						SHAPE(interiorElementBoundaries)[0],
						SHAPE(exteriorElementBoundaries)[0],
						SHAPE(dS)[1],
						SHAPE(dS)[2],
						SHAPE(elementResidual)[1],
						SHAPE(velocity)[2],
						IDATA(interiorElementBoundaries),
						IDATA(exteriorElementBoundaries),
						IDATA(elementBoundaryElements),
						IDATA(elementBoundaryLocalElementBoundaries),
						IDATA(exteriorElementBoundariesToSkip),
						DDATA(dS),
						DDATA(n),
						DDATA(elementResidual),
						DDATA(velocity),
						DDATA(conservationResidual));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cpostprocessingUpdateSelectedExteriorElementBoundaryFlux(PyObject* self, 
							 PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *skipflag_elementBoundaries,
    *flux,
    *w,
    *residual;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &skipflag_elementBoundaries,
                       &flux,
                       &w,
                       &residual))
    return NULL;
  updateSelectedExteriorElementBoundaryFlux(SHAPE(exteriorElementBoundaries)[0],
					    SHAPE(w)[1],
					    SHAPE(w)[2],
					    SHAPE(w)[3],
					    IDATA(exteriorElementBoundaries),
					    IDATA(elementBoundaryElements),
					    IDATA(elementBoundaryLocalElementBoundaries),
					    IDATA(skipflag_elementBoundaries),
					    DDATA(flux),
					    DDATA(w),
					    DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingUpdateAdvectiveVelocityPointEval(PyObject* self,
					  PyObject* args)
{
  PyObject *velocity,*advectiveVelocity;
  double updateCoef;
  int i,nPoints=1;
  if(!PyArg_ParseTuple(args,"dOO",
		       &updateCoef,
		       &advectiveVelocity,
		       &velocity))
    return NULL;
  for (i=0; i < ND(velocity)-1; i++)
    nPoints*=SHAPE(velocity)[i];
  postprocessAdvectiveVelocityPointEval(nPoints,
					SHAPE(velocity)[ND(velocity)-1],
					updateCoef,
					DDATA(advectiveVelocity),
					DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingUpdateDiffusiveVelocityPointEval(PyObject* self,
						PyObject* args)
{
  PyObject *velocity,*diffusionTensor,*grad_phi;
  double updateCoef;
  int i,nPoints=1;
  if(!PyArg_ParseTuple(args,"dOOO",
		       &updateCoef,
		       &diffusionTensor,
		       &grad_phi,
		       &velocity))
    return NULL;
  for (i=0; i < ND(velocity)-1; i++)
    nPoints*=SHAPE(velocity)[i];

  postprocessDiffusiveVelocityPointEval(nPoints,
					SHAPE(velocity)[ND(velocity)-1],
					updateCoef,
					DDATA(diffusionTensor),
					DDATA(grad_phi),
					DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingUpdateDiffusiveVelocityPointEval_sd(PyObject* self,
						   PyObject* args)
{
  PyObject *velocity,*diffusionTensor,*grad_phi,*rowptr,*colind;
  double updateCoef;
  int i,nPoints=1;
  if(!PyArg_ParseTuple(args,"dOOOOO",
		       &updateCoef,
		       &rowptr,
		       &colind,
		       &diffusionTensor,
		       &grad_phi,
		       &velocity))
    return NULL;
  for (i=0; i < ND(velocity)-1; i++)
    nPoints*=SHAPE(velocity)[i];

  postprocessDiffusiveVelocityPointEval_sd(nPoints,
					   SHAPE(velocity)[ND(velocity)-1],
					   updateCoef,
					   IDATA(rowptr),
					   IDATA(colind),
					   DDATA(diffusionTensor),
					   DDATA(grad_phi),
					   DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingCalculateElementResidualPWL(PyObject* self,
					   PyObject* args)
{
  PyObject *alpha,*elementResidual,*elementResidualPWL;
  if(!PyArg_ParseTuple(args,"OOO",
		       &alpha,
		       &elementResidual,
		       &elementResidualPWL))
    return NULL;
  calculateElementResidualPWL(SHAPE(elementResidual)[0],
			      SHAPE(elementResidual)[1],
			      SHAPE(elementResidualPWL)[1],
			      DDATA(alpha),
			      DDATA(elementResidual),
			      DDATA(elementResidualPWL));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingCopyElementBoundaryVelocityToParVec(PyObject* self,
						   PyObject* args)
{
  int ebN,k,I,nElementBoundaries,nQuadraturePoints,nSpace;
  PyObject *ebq_velocityO, *permutationsO, *ebq_v_par_localO;
  double *ebq_velocity,*ebq_v_par_local; 
  int *permutations;
  if(!PyArg_ParseTuple(args,"OOO",
		       &ebq_velocityO, 
		       &permutationsO, 
		       &ebq_v_par_localO))
    return NULL;
  ebq_velocity = DDATA(ebq_velocityO);
  ebq_v_par_local = DDATA(ebq_v_par_localO); 
  permutations = IDATA(permutationsO);
  nElementBoundaries=SHAPE(ebq_v_par_localO)[0];
  nQuadraturePoints=SHAPE(ebq_v_par_localO)[1];
  nSpace=SHAPE(ebq_v_par_localO)[2];
  for(ebN=0;ebN<nElementBoundaries;ebN++)
    for(k=0;k<nQuadraturePoints;k++)
      for(I=0;I<nSpace;I++)
	ebq_v_par_local[ebN*nQuadraturePoints*nSpace+permutations[ebN*nQuadraturePoints+k]*nSpace + I] = ebq_velocity[ebN*nQuadraturePoints*nSpace+k*nSpace + I];
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingAddAverageToParVec(PyObject* self,
				  PyObject* args)
{
  int ebN,k,I,nElementBoundaries,nQuadraturePoints,nSpace;
  PyObject *ebq_velocityAverageO, *permutationsO, *ebq_v_par_localO;
  double *ebq_velocityAverage,*ebq_v_par_local; 
  int *permutations;
  if(!PyArg_ParseTuple(args,"OOO",
		       &ebq_velocityAverageO, 
		       &permutationsO, 
		       &ebq_v_par_localO))
    return NULL;
  ebq_velocityAverage = DDATA(ebq_velocityAverageO);
  ebq_v_par_local = DDATA(ebq_v_par_localO); 
  permutations = IDATA(permutationsO);
  nElementBoundaries=SHAPE(ebq_v_par_localO)[0];
  nQuadraturePoints=SHAPE(ebq_v_par_localO)[1];
  nSpace=SHAPE(ebq_v_par_localO)[2];
  for(ebN=0;ebN<nElementBoundaries;ebN++)
    for(k=0;k<nQuadraturePoints;k++)
      for(I=0;I<nSpace;I++)
	ebq_v_par_local[ebN*nQuadraturePoints*nSpace+permutations[ebN*nQuadraturePoints+k]*nSpace + I] += ebq_velocityAverage[ebN*nQuadraturePoints*nSpace+k*nSpace + I];
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cpostprocessingCopyParVecToElementBoundaryVelocity(PyObject* self,
						   PyObject* args)
{
  int ebN,k,I,nElementBoundaries,nQuadraturePoints,nSpace;
  PyObject *ebq_velocityO, *permutationsO, *ebq_v_par_localO;
  double *ebq_velocity,*ebq_v_par_local; 
  int *permutations;
  if(!PyArg_ParseTuple(args,"OOO",
		       &ebq_velocityO, 
		       &permutationsO, 
		       &ebq_v_par_localO))
    return NULL;
  ebq_velocity = DDATA(ebq_velocityO);
  ebq_v_par_local = DDATA(ebq_v_par_localO); 
  permutations = IDATA(permutationsO);
  nElementBoundaries=SHAPE(ebq_v_par_localO)[0];
  nQuadraturePoints=SHAPE(ebq_v_par_localO)[1];
  nSpace=SHAPE(ebq_v_par_localO)[2];
  for(ebN=0;ebN<nElementBoundaries;ebN++)
    for(k=0;k<nQuadraturePoints;k++)
      for(I=0;I<nSpace;I++)
	ebq_velocity[ebN*nQuadraturePoints*nSpace+k*nSpace + I] = ebq_v_par_local[ebN*nQuadraturePoints*nSpace+permutations[ebN*nQuadraturePoints+k]*nSpace + I];
  Py_INCREF(Py_None);
  return Py_None;
}


static PyMethodDef cpostprocessingMethods[] = {
  { "postProcessRT0velocityFromP1nc", 
    cpostprocessingPostProcessRT0velocityFromP1nc,
    METH_VARARGS, 
    "convert P1nc flux to RT0 flux locally allow for full mass matrix, etc"},
  { "postProcessRT0velocityFromP1nc_sd", 
    cpostprocessingPostProcessRT0velocityFromP1nc_sd,
    METH_VARARGS, 
    "convert P1nc flux to RT0 flux locally allow for full mass matrix, etc"},
  { "getElementRT0velocityValues", 
    cpostprocessingGetElementRT0velocityValues,
    METH_VARARGS, 
    "calculate P1nc flux values at points defined on elements"},
  { "updateRT0velocityWithAveragedPotentialP1nc", 
    cpostprocessingUpdateRT0velocityWithAveragedPotentialP1nc,
    METH_VARARGS, 
    "incorporate additional constant flux terms for multiple potentials"},
  { "updateRT0velocityWithAveragedPotentialP1nc_sd", 
    cpostprocessingUpdateRT0velocityWithAveragedPotentialP1nc_sd,
    METH_VARARGS, 
    "incorporate additional constant flux terms for multiple potentials"},
  { "getElementBoundaryRT0velocityValues", 
    cpostprocessingGetElementBoundaryRT0velocityValues,
    METH_VARARGS, 
    "calculate P1nc flux values at points defined on element boundaries"},
  { "getGlobalElementBoundaryRT0velocityValues", 
    cpostprocessingGetGlobalElementBoundaryRT0velocityValues,
    METH_VARARGS, 
    "calculate P1nc flux values at points defined uniquely on global exterior elementBoundaries"},
  { "getGlobalExteriorElementBoundaryRT0velocityValues", 
    cpostprocessingGetGlobalExteriorElementBoundaryRT0velocityValues,
    METH_VARARGS, 
    "calculate P1nc flux values at points defined uniquely on elementBoundaries"},
  { "getElementBoundaryRT0velocityValuesFluxRep", 
    cpostprocessingGetElementBoundaryRT0velocityValuesFluxRep,
    METH_VARARGS, 
    "calculate P1nc flux values at points defined on element boundaries"},
  { "getGlobalElementBoundaryRT0velocityValuesFluxRep", 
    cpostprocessingGetGlobalElementBoundaryRT0velocityValuesFluxRep,
    METH_VARARGS, 
    "calculate P1nc flux values at points defined uniquely on elementBoundaries"},
  { "getGlobalExteriorElementBoundaryRT0velocityValuesFluxRep", 
    cpostprocessingGetGlobalExteriorElementBoundaryRT0velocityValuesFluxRep,
    METH_VARARGS, 
    "calculate P1nc flux values at points defined uniquely on exterior elementBoundaries"},
  { "getRT0velocityValuesFluxRep_arbitraryElementMembership", 
    cpostprocessingGetRT0velocityValuesFluxRep_arbitraryElementMembership,
    METH_VARARGS, 
    "calculate RT0 velocity at points {x} located given index of element locations"},
  {"fluxCorrectionVelocityUpdate",
   cpostprocessingFluxCorrectionVelocityUpdate,
   METH_VARARGS,
   "update velocity for global postprocessing (pwc, Sun-Wheeler)"},
  {"computeFluxCorrectionPWC",
   cpostprocessingComputeFluxCorrectionPWC,
   METH_VARARGS,
   "calculate flux correction on element boundary using pwc element variables"},
  {"sunWheelerGSsweep",
   cpostprocessingSunWheelerGSsweep,
   METH_VARARGS,
   "calculate one sweep for Sun-Wheeler iteration"},
  {"calculateConservationResidualGlobalBoundaries",
   cpostprocessingCalculateConservationResidualGlobalBoundaries,
   METH_VARARGS,
   "calculate the conservation residuals using global boundary info"},
  {"calculateConservationResidualPWL",
   cpostprocessingCalculateConservationResidualPWL,
   METH_VARARGS,
   "calculate the conservation residuals for star shaped subdomains"},
  {"calculateConservationJacobianPWL",
   cpostprocessingCalculateConservationJacobianPWL,
   METH_VARARGS,
   "calculate the (LU-factorized) conservation jacobians for  star shaped subdomains"},
  {"calculateConservationFluxPWL",
   cpostprocessingCalculateConservationFluxPWL,
   METH_VARARGS,
   "solve the mass conservation system and calculate the conservative flux and residual"},
  {"calculateConservationFluxPWL_noNeumannFix",
   cpostprocessingCalculateConservationFluxPWL_noNeumannFix,
   METH_VARARGS,
   "solve the mass conservation system and calculate the conservative flux and residual -- does not include manual fix for pure Neumann problems"},
  {"calculateConservationResidualPWL_opt",
   cpostprocessingCalculateConservationResidualPWL_opt,
   METH_VARARGS,
   "calculate the conservation residuals for star shaped subdomains"},
  {"calculateConservationResidualPWL_primative",
   cpostprocessingCalculateConservationResidualPWL_primative,
   METH_VARARGS,
   "calculate the conservation residuals for star shaped subdomains"},
  {"calculateConservationJacobianPWL_opt",
   cpostprocessingCalculateConservationJacobianPWL_opt,
   METH_VARARGS,
   "calculate the (LU-factorized) conservation jacobians for  star shaped subdomains"},
  {"calculateConservationFluxPWL_opt",
   cpostprocessingCalculateConservationFluxPWL_opt,
   METH_VARARGS,
   "solve the mass conservation system and calculate the conservative flux and residual"},
  {"calculateConservationResidualPWL_interiorBoundaries",
   cpostprocessingCalculateConservationResidualPWL_interiorBoundaries,
   METH_VARARGS,
   "calculate the conservation residuals for star shaped subdomains -- allows internal Neumann boundaries"},
  {"calculateConservationJacobianPWL_interiorBoundaries",
   cpostprocessingCalculateConservationJacobianPWL_interiorBoundaries,
   METH_VARARGS,
   "calculate the (LU-factorized) conservation jacobians for star shaped subdomains -- allows internal Neumann boundaries"},
  { "postProcessRT0potentialFromP1nc", 
    cpostprocessingPostProcessRT0potentialFromP1nc,
    METH_VARARGS, 
    "convert P1nc potential to RT0 potential locally"},
  { "postProcessRT0potentialFromP1nc_sd", 
    cpostprocessingPostProcessRT0potentialFromP1nc_sd,
    METH_VARARGS, 
    "convert P1nc potential to RT0 potential locally"},
  { "projectElementBoundaryVelocityToRT0fluxRep", 
    cpostprocessingProjectElementBoundaryVelocityToRT0fluxRep,
    METH_VARARGS, 
    "project velocity defined on element boundary to local RT0 basis in 'flux' rep"},
  { "projectElementBoundaryFluxToRT0fluxRep", 
    cpostprocessingProjectElementBoundaryFluxToRT0fluxRep,
    METH_VARARGS, 
    "project normal flux defined on element boundary to local RT0 basis in 'flux' rep"},
  { "getElementRT0velocityValuesFluxRep", 
    cpostprocessingGetElementRT0velocityValuesFluxRep,
    METH_VARARGS, 
    "get velocity at element quad points assuming defined using local RT0 basis in 'flux' rep"},
  { "buildLocalBDM1projectionMatrices", 
    cpostprocessingBuildLocalBDM1projectionMatrices,
    METH_VARARGS, 
    "build local projections to BDM1 assuming standard basis for P^1(E)"},
  { "buildLocalBDM2projectionMatrices", 
    cpostprocessingBuildLocalBDM2projectionMatrices,
    METH_VARARGS, 
    "build local projections to BDM2 assuming standard basis for P^2(E)"},
  { "factorLocalBDM1projectionMatrices", 
    cpostprocessingFactorLocalBDM1projectionMatrices,
    METH_VARARGS, 
    "compute LU factorization for local BDM projection matrices"},
  { "factorLocalBDM2projectionMatrices", 
    cpostprocessingFactorLocalBDM2projectionMatrices,
    METH_VARARGS, 
    "compute LU factorization for local BDM projection matrices"},
  { "solveLocalBDM1projection", 
    cpostprocessingSolveLocalBDM1projection,
    METH_VARARGS, 
    "solve for local [P^1(E)]^d dofs using BDM1 projection assuming system already factored" },
  { "solveLocalBDM2projection", 
    cpostprocessingSolveLocalBDM2projection,
    METH_VARARGS, 
    "solve for local [P^1(E)]^d dofs using BDM1 projection assuming system already factored" },
  { "buildBDM2rhs", 
    cpostprocessingBuildBDM2rhs,
    METH_VARARGS, 
    "solve for local [P^1(E)]^d dofs using BDM2 projection assuming system already factored" },
  { "solveLocalBDM1projectionFromFlux", 
    cpostprocessingSolveLocalBDM1projectionFromFlux,
    METH_VARARGS, 
    "solve for local [P^1(E)]^d dofs using BDM1 projection assuming system already factored" },
  { "getElementBDM1velocityValuesLagrangeRep", 
    cpostprocessingGetElementBDM1velocityValuesLagrangeRep,
    METH_VARARGS, 
    "get velocity at quadrature points assuming have local dofs for local P^1(E) and std Lagr. basis" },
  { "getElementBDM2velocityValuesLagrangeRep", 
    cpostprocessingGetElementBDM2velocityValuesLagrangeRep,
    METH_VARARGS, 
    "get velocity at quadrature points assuming have local dofs for local P^1(E) and std Lagr. basis" },
  { "getElementLDGvelocityValuesLagrangeRep", 
    cpostprocessingGetElementLDGvelocityValuesLagrangeRep,
    METH_VARARGS, 
    "get velocity at quadrature points assuming have local dofs for local P^k(E) where k=1,2 and std Lagr. basis" },
  { "getGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep", 
    cpostprocessingGetGlobalExteriorElementBoundaryBDM1velocityValuesLagrangeRep,
    METH_VARARGS, 
    "get velocity at exterior element boundary quadrature points assuming have local dofs for local P^1(E) and std Lagr. basis" },
 { "getGlobalElementBoundaryBDM1velocityValuesLagrangeRep", 
    cpostprocessingGetGlobalElementBoundaryBDM1velocityValuesLagrangeRep,
    METH_VARARGS, 
    "get velocity at exterior element boundary quadrature points (but stored in the global element boundary quadrature array) assuming have local dofs for local P^1(E) and std Lagr. basis" },
  { "getElementBoundaryBDM1velocityValuesLagrangeRep", 
    cpostprocessingGetElementBoundaryBDM1velocityValuesLagrangeRep,
    METH_VARARGS, 
    "get velocity at element boundary quadrature points assuming have local dofs for local P^1(E) and std Lagr. basis" },
  { "updateSelectedExteriorElementBoundaryFlux", 
    cpostprocessingUpdateSelectedExteriorElementBoundaryFlux,
    METH_VARARGS, 
    "update element residual using flux from flagged element boundaries only" },
  { "updateAdvectiveVelocityPointEval", 
    cpostprocessingUpdateAdvectiveVelocityPointEval,
    METH_VARARGS, 
    "evaluate advective velocity using pointwise definition" },
  { "updateDiffusiveVelocityPointEval", 
    cpostprocessingUpdateDiffusiveVelocityPointEval,
    METH_VARARGS, 
    "evaluate advective velocity using pointwise definition" },
  { "updateDiffusiveVelocityPointEval_sd", 
    cpostprocessingUpdateDiffusiveVelocityPointEval_sd,
    METH_VARARGS, 
    "evaluate advective velocity using pointwise definition" },
  { "subdomain_U_copy_local2global", 
    cpostprocessing_subdomain_U_copy_local2global,
    METH_VARARGS, 
    "copy over the corrections from node star factors to parallel storage" },
  { "subdomain_U_copy_global2local", 
    cpostprocessing_subdomain_U_copy_global2local,
    METH_VARARGS, 
    "copy over the corrections from parallel storage to node star factors" },
  { "calculateElementResidualPWL", 
    cpostprocessingCalculateElementResidualPWL,
    METH_VARARGS, 
    "calculate the element residual for the linear nodal (Lagrange) basis from the element residuals in a different basis" },
  { "copyElementBoundaryVelocityToParVec", 
    cpostprocessingCopyElementBoundaryVelocityToParVec,
    METH_VARARGS, 
    "Copy the velocity into the local parallel (ghosted) storage" },
  { "addAverageToParVec",
    cpostprocessingAddAverageToParVec,
    METH_VARARGS, 
    "Updated the local parallel (ghosted) storage with the average" },
  { "copyParVecToElementBoundaryVelocity",
    cpostprocessingCopyParVecToElementBoundaryVelocity,
    METH_VARARGS, 
    "Copy the parallel storage into the velocity" },
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcpostprocessing(void)
{
  PyObject *m,*d;
  if (PyType_Ready(&NodeStarFactorType) < 0)
    return;
  m = Py_InitModule3("cpostprocessing", cpostprocessingMethods,"postprocessing module");

  Py_INCREF(&NodeStarFactorType);
  PyModule_AddObject(m, "NodeStarFactor", (PyObject *)&NodeStarFactorType);
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
