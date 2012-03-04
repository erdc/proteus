#include "Python.h"
#include "numpy/arrayobject.h"
#include "ellam.h"
#include "superluWrappersModule.h"
#include <vector>
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

static PyObject* cellam_updateOldMass_weak(PyObject* self,
					   PyObject* args)
{
  int nSpace,            
    nDOF_test_element, 
    nElements_global,  
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nQuadraturePoints_element;
  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *elementBoundaryOuterNormalsArray,
    *dV,
    *x_track,
    *t_track,
    *element_track,
    *flag_track,
    *u_l2g, 
    *q_m_last,
    *elementResidual_u; 

  /*return value*/
  double totalOldMass = 0.0;

  if (!PyArg_ParseTuple(args,
                        "iiiiiiiOOOOOOOOOOOO",/*iiO",*/
			&nSpace,
			&nDOF_test_element,
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
			&u_l2g, 
			&q_m_last,
			&elementResidual_u))

    return NULL;
  
  totalOldMass = updateOldMass_weak(nSpace,
				    nDOF_test_element,
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
				    IDATA(u_l2g), 
				    DDATA(q_m_last),
				    DDATA(elementResidual_u));
  
  return Py_BuildValue("d",totalOldMass);
}
static PyObject* cellam_updateOldMass_weak_arbitraryQuadrature(PyObject* self,
							       PyObject* args)
{
  int nSpace,            
    nDOF_test_element, 
    nElements_global,  
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nQuadraturePoints_track;
  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *elementBoundaryOuterNormalsArray,
    *dV_track,
    *x_track,
    *t_track,
    *element_track,
    *flag_track,
    *u_l2g, 
    *q_m_track,
    *elementResidual_u; 

  if (!PyArg_ParseTuple(args,
                        "iiiiiiiOOOOOOOOOOOO",/*iiO",*/
			&nSpace,
			&nDOF_test_element,
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nQuadraturePoints_track,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&elementBoundaryOuterNormalsArray,
			&dV_track,
			&x_track,
			&t_track,
			&element_track,
			&flag_track,
			&u_l2g, 
			&q_m_track,
			&elementResidual_u))

    return NULL;
  
  updateOldMass_weak_arbitraryQuadrature(nSpace,
					 nDOF_test_element,
					 nElements_global,
					 nNodes_global,
					 nNodes_element,
					 nElementBoundaries_element,
					 nQuadraturePoints_track,
					 DDATA(nodeArray),
					 IDATA(elementNodesArray),
					 IDATA(elementNeighborsArray),
					 DDATA(elementBoundaryOuterNormalsArray),
					 DDATA(dV_track),
					 DDATA(x_track),
					 DDATA(t_track),
					 IDATA(element_track),
					 IDATA(flag_track),
					 IDATA(u_l2g), 
					 DDATA(q_m_track),
					 DDATA(elementResidual_u));
  
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_updateNewMass_weak(PyObject* self,
					   PyObject* args)
{
  int nSpace,            
    nDOF_test_element, 
    nElements_global,  
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
    nQuadraturePoints_element;
  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *elementBoundaryOuterNormalsArray,
    *dV,
    *x,
    *u_l2g, 
    *q_m,
    *elementResidual_u; 

  /*return value*/
  double totalNewMass = 0.0;

  if (!PyArg_ParseTuple(args,
                        "iiiiiiiOOOOOOOOO",/*iiO",*/
			&nSpace,
			&nDOF_test_element,
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
			&x,
			&u_l2g, 
			&q_m,
			&elementResidual_u))

    return NULL;
  
  totalNewMass = updateNewMass_weak(nSpace,
				    nDOF_test_element,
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
				    DDATA(x),
				    IDATA(u_l2g), 
				    DDATA(q_m),
				    DDATA(elementResidual_u));
  
  return Py_BuildValue("d",totalNewMass);
}

static PyObject* cellam_evaluateSolutionAtTrackedPoints(PyObject* self,
							PyObject* args)
{
  int nSpace,            
    nDOF_trial_element, 
    nPoints_tracked,
    nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element;
  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *elementBoundaryOuterNormalsArray,
    *x_track,
    *t_track,
    *element_track,
    *flag_track,
    *u_l2g, 
    *u_dof,
    *u_x_track;

  if (!PyArg_ParseTuple(args,
                        "iiiiiiiOOOOOOOOOOO",
			&nSpace,
			&nDOF_trial_element,
			&nPoints_tracked,
			&nElements_global,
			&nNodes_global,
			&nNodes_element,
			&nElementBoundaries_element,
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&elementBoundaryOuterNormalsArray,
			&x_track,
			&t_track,
			&element_track,
			&flag_track,
			&u_l2g, 
			&u_dof,
			&u_x_track))
    return NULL;
  
  evaluateSolutionAtTrackedPoints(nSpace,
				  nDOF_trial_element,
				  nPoints_tracked,
				  nElements_global,
				  nNodes_global,
				  nNodes_element,
				  nElementBoundaries_element,
				  DDATA(nodeArray),
				  IDATA(elementNodesArray),
				  IDATA(elementNeighborsArray),
				  DDATA(elementBoundaryOuterNormalsArray),
				  DDATA(x_track),
				  DDATA(t_track),
				  IDATA(element_track),
				  IDATA(flag_track),
				  IDATA(u_l2g), 
				  DDATA(u_dof),
				  DDATA(u_x_track)); 
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cellam_updateExteriorOutflowBoundaryFlux(PyObject* self,
							  PyObject* args)
{
  double dtnp1;
  int nSpace,
    nDOF_test_element,
    nQuadraturePoints_elementBoundary,
    nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *ebqe_velocity_ext,
    *ebqe_n_ext,
    *ebqe_outflow_flux_last,
    *u_test_dS_ext,
    *ebqe_u,
    *u_l2g,
    *ebqe_outflow_flux;
  PyObject* q_elementResidual_u;
  /*return value */
  double totalOutflowFlux = 0.0;

  if (!PyArg_ParseTuple(args,
                        "diiiiOOOOOOOOOOO",
			&dtnp1,
			&nSpace,
			&nDOF_test_element,
			&nQuadraturePoints_elementBoundary,
			&nExteriorElementBoundaries_global,
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
			&ebqe_velocity_ext,
			&ebqe_n_ext,
			&ebqe_outflow_flux_last,
			&u_test_dS_ext,
			&ebqe_u,
			&u_l2g,
			&ebqe_outflow_flux,
			&q_elementResidual_u)) 
			
    return NULL;

  totalOutflowFlux = updateExteriorOutflowBoundaryFlux(dtnp1,
						       nSpace,
						       nDOF_test_element,
						       nQuadraturePoints_elementBoundary,
						       nExteriorElementBoundaries_global,
						       IDATA(exteriorElementBoundariesArray),
						       IDATA(elementBoundaryElementsArray),
						       IDATA(elementBoundaryLocalElementBoundariesArray),
						       DDATA(ebqe_velocity_ext),
						       DDATA(ebqe_n_ext),
						       DDATA(ebqe_outflow_flux_last),
						       DDATA(u_test_dS_ext),
						       DDATA(ebqe_u),
						       IDATA(u_l2g),
						       DDATA(ebqe_outflow_flux),
						       DDATA(q_elementResidual_u)); 

  return Py_BuildValue("d",totalOutflowFlux);

} 

static PyObject* cellam_updateExteriorOutflowBoundaryFluxInGlobalResidual(PyObject* self,
									  PyObject* args)
{
  double dtnp1;
  int nSpace,
    nDOF_test_element,
    nQuadraturePoints_elementBoundary,
    nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *ebqe_velocity_ext,
    *ebqe_n_ext,
    *ebqe_outflow_flux_last,
    *u_test_dS_ext,
    *ebqe_u,
    *u_l2g,
    *ebqe_outflow_flux;
  int offset_u,stride_u;
  PyObject* q_elementResidual_u,
    *globalResidual;

  /*return value*/
  double totalOutflowFlux = 0.0;
  if (!PyArg_ParseTuple(args,
                        "diiiiOOOOOOOOOOiiOO",
			&dtnp1,
			&nSpace,
			&nDOF_test_element,
			&nQuadraturePoints_elementBoundary,
			&nExteriorElementBoundaries_global,
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
			&ebqe_velocity_ext,
			&ebqe_n_ext,
			&ebqe_outflow_flux_last,
			&u_test_dS_ext,
			&ebqe_u,
			&u_l2g,
			&ebqe_outflow_flux,
			&offset_u, &stride_u, 
			&q_elementResidual_u, 
			&globalResidual))
    return NULL;

  totalOutflowFlux = updateExteriorOutflowBoundaryFluxInGlobalResidual(dtnp1,
								       nSpace,
								       nDOF_test_element,
								       nQuadraturePoints_elementBoundary,
								       nExteriorElementBoundaries_global,
								       IDATA(exteriorElementBoundariesArray),
								       IDATA(elementBoundaryElementsArray),
								       IDATA(elementBoundaryLocalElementBoundariesArray),
								       DDATA(ebqe_velocity_ext),
								       DDATA(ebqe_n_ext),
								       DDATA(ebqe_outflow_flux_last),
								       DDATA(u_test_dS_ext),
								       DDATA(ebqe_u),
								       IDATA(u_l2g),
								       DDATA(ebqe_outflow_flux),
								       offset_u,stride_u, 
								       DDATA(q_elementResidual_u), 
								       DDATA(globalResidual));


  return Py_BuildValue("d",totalOutflowFlux);

} 


static PyObject* cellam_updateExteriorOutflowBoundaryFluxGlobalJacobian(PyObject* self,
									PyObject* args)
{
  double dtnp1;
  int nSpace,
    nDOF_test_element,
    nDOF_trial_element,
    nQuadraturePoints_elementBoundary,
    nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *ebqe_velocity_ext,
    *ebqe_n_ext,
    *ebqe_outflow_flux_last,
    *u_test_dS_ext,
    *ebqe_u,
    *u_trial_ext,
    *csrRowIndeces_u_u, *csrColumnOffsets_u_u,
    *csrColumnOffsets_eb_u_u,
    *globalJacobian;


  if (!PyArg_ParseTuple(args,
                        "diiiiiOOOOOOOOOOOOO",
			&dtnp1,
			&nSpace,
			&nDOF_test_element,
			&nDOF_trial_element,
			&nQuadraturePoints_elementBoundary,
			&nExteriorElementBoundaries_global,
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
			&ebqe_velocity_ext,
			&ebqe_n_ext,
			&ebqe_outflow_flux_last,
			&u_test_dS_ext,
			&ebqe_u,
			&u_trial_ext,
			&csrRowIndeces_u_u,
			&csrColumnOffsets_u_u,
			&csrColumnOffsets_eb_u_u,
			&globalJacobian))
    return NULL;

  updateExteriorOutflowBoundaryFluxGlobalJacobian(dtnp1,
						  nSpace,
						  nDOF_test_element,
						  nDOF_trial_element,
						  nQuadraturePoints_elementBoundary,
						  nExteriorElementBoundaries_global,
						  IDATA(exteriorElementBoundariesArray),
						  IDATA(elementBoundaryElementsArray),
						  IDATA(elementBoundaryLocalElementBoundariesArray),
						  DDATA(ebqe_velocity_ext),
						  DDATA(ebqe_n_ext),
						  DDATA(ebqe_outflow_flux_last),
						  DDATA(u_test_dS_ext),
						  DDATA(ebqe_u),
						  DDATA(u_trial_ext),
						  IDATA(csrRowIndeces_u_u),
						  IDATA(csrColumnOffsets_u_u),
						  IDATA(csrColumnOffsets_eb_u_u),
						  CSRVAL(globalJacobian));
  
  Py_INCREF(Py_None); 
  return Py_None;

} 

static PyObject* cellam_updateExteriorOutflowBoundaryFluxJacobian(PyObject* self,
								  PyObject* args)
{
  double dtnp1;
  int nSpace,
    nDOF_test_element,
    nDOF_trial_element,
    nQuadraturePoints_elementBoundary,
    nExteriorElementBoundaries_global;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *ebqe_velocity_ext,
    *ebqe_n_ext,
    *ebqe_outflow_flux_last,
    *u_test_dS_ext,
    *ebqe_u,
    *u_trial_ext,
    *fluxJacobian;


  if (!PyArg_ParseTuple(args,
                        "diiiiiOOOOOOOOOO",
			&dtnp1,
			&nSpace,
			&nDOF_test_element,
			&nDOF_trial_element,
			&nQuadraturePoints_elementBoundary,
			&nExteriorElementBoundaries_global,
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
			&ebqe_velocity_ext,
			&ebqe_n_ext,
			&ebqe_outflow_flux_last,
			&u_test_dS_ext,
			&ebqe_u,
			&u_trial_ext,
			&fluxJacobian))
    return NULL;

  updateExteriorOutflowBoundaryFluxJacobian(dtnp1,
					    nSpace,
					    nDOF_test_element,
					    nDOF_trial_element,
					    nQuadraturePoints_elementBoundary,
					    nExteriorElementBoundaries_global,
					    IDATA(exteriorElementBoundariesArray),
					    IDATA(elementBoundaryElementsArray),
					    IDATA(elementBoundaryLocalElementBoundariesArray),
					    DDATA(ebqe_velocity_ext),
					    DDATA(ebqe_n_ext),
					    DDATA(ebqe_outflow_flux_last),
					    DDATA(u_test_dS_ext),
					    DDATA(ebqe_u),
					    DDATA(u_trial_ext),
					    DDATA(fluxJacobian));
  
  Py_INCREF(Py_None); 
  return Py_None;

} 

static PyObject* cellam_markInflowBoundaryPoints(PyObject* self,
						 PyObject* args)
{
  int nSpace;

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
                        "idddiiOOOOOOOOOOO",
			&nSpace,
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

  markInflowBoundaryPoints(nSpace,
			   tn,
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

static PyObject* cellam_accumulateInflowFluxInGlobalResidual(PyObject* self,
							     PyObject* args)
{
  int nSpace,
    nDOF_test_element;

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

  /*return value is total accumulated inflow flux*/
  double totalInflowFlux;
  if (!PyArg_ParseTuple(args,
			"iiiiiiiiOOOOOOOddOOOOOOOOOOiiOOO",
			&nSpace,
			&nDOF_test_element,
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

  totalInflowFlux = accumulateInflowFluxInGlobalResidual(nSpace,
							 nDOF_test_element,
							 nElements_global,
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


  return Py_BuildValue("d",totalInflowFlux);
}

static PyObject* cellam_accumulateInflowFlux(PyObject* self,
					     PyObject* args)
{
  int nSpace,
    nDOF_test_element;

  double tp,timeWeight;

  int nElements_global,
    nNodes_global,
    nNodes_element,
    nElementBoundaries_element,
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
    *sdInfo_u_rowptr,
    *sdInfo_u_colind,
    *isFluxBoundary_u,
    *ebqe_bc_flux_u_ext;

  /*return value is total accumulated inflow flux*/
  double totalInflowFlux;

  if (!PyArg_ParseTuple(args,
			"iiiiiiiiOOOOOOOddOOOOOOOOOOOO",
			&nSpace,
			&nDOF_test_element,
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
			&isFluxBoundary_u,
			&ebqe_bc_flux_u_ext))
    return NULL;

  totalInflowFlux = accumulateInflowFlux(nSpace,
					 nDOF_test_element,
					 nElements_global,
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
					 IDATA(isFluxBoundary_u),
					 DDATA(ebqe_bc_flux_u_ext));

  return Py_BuildValue("d",totalInflowFlux);
 
}


static PyObject* cellam_tagNegligibleIntegrationPoints(PyObject* self,
						       PyObject* args)
{
  int nPoints;
  double zeroTol;
  PyObject *x,*u,*flag_track;

  if (!PyArg_ParseTuple(args,"idOOO",
			&nPoints,
			&zeroTol,
			&x,
			&u,
			&flag_track))
    return NULL;

  tagNegligibleIntegrationPoints(nPoints,
				 zeroTol,
				 DDATA(x),
				 DDATA(u),
				 IDATA(flag_track));

  
  Py_INCREF(Py_None); 
  return Py_None;
 
}
static PyObject* cellam_calculateSlumpedMassApproximation1d(PyObject* self,
							    PyObject* args)
{
  PyObject *u_l2g,
    *elementNeighborsArray,
    *u_dof,*u_dof_limit,
    *dm,*w,*v,
    *dV,*rhs,
    *theta,
    *slumpedMassMatrixCorrection,
    *elementResidual; 

  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOOO",
			&u_l2g,
			&elementNeighborsArray,
			&u_dof,
			&u_dof_limit,
			&dm,
			&w,
			&v,
			&dV,
			&rhs,
			&elementResidual,
			&theta,
			&slumpedMassMatrixCorrection))
    return NULL;
  
  calculateSlumpedMassApproximation1d(SHAPE(u_l2g)[0],
				      SHAPE(elementNeighborsArray)[1],
				      SHAPE(dm)[1],
				      SHAPE(v)[2],
				      SHAPE(w)[2],
				      IDATA(u_l2g),
				      IDATA(elementNeighborsArray),
				      DDATA(u_dof),
				      DDATA(u_dof_limit),
				      DDATA(dm),
				      DDATA(w),
				      DDATA(v),
				      DDATA(dV),
				      DDATA(rhs),
				      DDATA(elementResidual),
				      DDATA(theta),
				      DDATA(slumpedMassMatrixCorrection));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_calculateSlumpedMassApproximation1d_local(PyObject* self,
								  PyObject* args)
{
  PyObject *u_l2g,
    *elementNeighborsArray,
    *u_dof,*u_dof_limit,
    *dm,*w,*v,
    *dV,*rhs,
    *theta,
    *slumpedMassMatrixCorrection,
    *elementResidual; 

  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOOO",
			&u_l2g,
			&elementNeighborsArray,
			&u_dof,
			&u_dof_limit,
			&dm,
			&w,
			&v,
			&dV,
			&rhs,
			&elementResidual,
			&theta,
			&slumpedMassMatrixCorrection))
    return NULL;
  
  calculateSlumpedMassApproximation1d_local(SHAPE(u_l2g)[0],
					    SHAPE(elementNeighborsArray)[1],
					    SHAPE(dm)[1],
					    SHAPE(v)[2],
					    SHAPE(w)[2],
					    IDATA(u_l2g),
					    IDATA(elementNeighborsArray),
					    DDATA(u_dof),
					    DDATA(u_dof_limit),
					    DDATA(dm),
					    DDATA(w),
					    DDATA(v),
					    DDATA(dV),
					    DDATA(rhs),
					    DDATA(elementResidual),
					    DDATA(theta),
					    DDATA(slumpedMassMatrixCorrection));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_calculateSlumpedMassApproximation2d(PyObject* self,
							    PyObject* args)
{
  PyObject *u_l2g,
    *elementNeighborsArray,
    *u_dof,*u_dof_limit,
    *dm,*w,*v,
    *dV,*rhs,
    *theta,
    *slumpedMassMatrixCorrection,
    *elementResidual; 
  double adjustFactor=1.0;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOOO|d",
			&u_l2g,
			&elementNeighborsArray,
			&u_dof,
			&u_dof_limit,
			&dm,
			&w,
			&v,
			&dV,
			&rhs,
			&elementResidual,
			&theta,
			&slumpedMassMatrixCorrection,
			&adjustFactor))
    return NULL;
  
  calculateSlumpedMassApproximation2d(SHAPE(u_l2g)[0],
				      SHAPE(elementNeighborsArray)[1],
				      SHAPE(dm)[1],
				      SHAPE(v)[2],
				      SHAPE(w)[2],
				      adjustFactor,
				      IDATA(u_l2g),
				      IDATA(elementNeighborsArray),
				      DDATA(u_dof),
				      DDATA(u_dof_limit),
				      DDATA(dm),
				      DDATA(w),
				      DDATA(v),
				      DDATA(dV),
				      DDATA(rhs),
				      DDATA(elementResidual),
				      DDATA(theta),
				      DDATA(slumpedMassMatrixCorrection));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_calculateSlumpedMassApproximation2d_upwind(PyObject* self,
								   PyObject* args)
{
  PyObject *nodeArray,
    *elementNodesArray,
    *elementNeighborsArray,
    *nodeStarOffsets,
    *nodeStarArray,
    *elementBoundaryLocalOuterNormalsArray,
    *u_l2g,*u_dof,*u_dof_limit,
    *dm,*df,*w,*v,
    *dV,*rhs,
    *theta,
    *slumpedMassMatrixCorrection,
    *elementResidual; 

  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOOOOOOOOO",
			&nodeArray,
			&elementNodesArray,
			&elementNeighborsArray,
			&nodeStarOffsets,
			&nodeStarArray,
			&elementBoundaryLocalOuterNormalsArray,
			&u_l2g,
			&u_dof,
			&u_dof_limit,
			&dm,
			&df,
			&w,
			&v,
			&dV,
			&rhs,
			&elementResidual,
			&theta,
			&slumpedMassMatrixCorrection))
    return NULL;
  
  calculateSlumpedMassApproximation2d_upwind(SHAPE(u_l2g)[0],
					     SHAPE(nodeArray)[0],
					     SHAPE(elementNeighborsArray)[1],
					     SHAPE(elementNodesArray)[1],
					     SHAPE(dm)[1],
					     SHAPE(v)[2],
					     SHAPE(w)[2],
					     DDATA(nodeArray),
					     IDATA(elementNodesArray),
					     IDATA(elementNeighborsArray),
					     IDATA(nodeStarOffsets),
					     IDATA(nodeStarArray),
					     DDATA(elementBoundaryLocalOuterNormalsArray),
					     IDATA(u_l2g),
					     DDATA(u_dof),
					     DDATA(u_dof_limit),
					     DDATA(dm),
					     DDATA(df),
					     DDATA(w),
					     DDATA(v),
					     DDATA(dV),
					     DDATA(rhs),
					     DDATA(elementResidual),
					     DDATA(theta),
					     DDATA(slumpedMassMatrixCorrection));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_updateElementJacobianWithSlumpedMassApproximation(PyObject* self,
									  PyObject* args)
{
  PyObject *theta,
    *elementJacobian; 

  if (!PyArg_ParseTuple(args,
                        "OO",
			&theta,
			&elementJacobian))
    return NULL;
  
  updateElementJacobianWithSlumpedMassApproximation(SHAPE(elementJacobian)[0],
						    SHAPE(elementJacobian)[2],
						    SHAPE(elementJacobian)[1],
						    DDATA(theta),
						    DDATA(elementJacobian));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cellam_calculateBerzinsSlumpedMassApproximation1d(PyObject* self,
								  PyObject* args)
{
  PyObject *u_l2g,
    *elementNeighborsArray,
    *u_dof,*u_dof_limit,
    *dm,*w,*v,
    *dV,*rhs,
    *slumpedMassMatrixCorrection,
    *elementResidual; 

  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOO",
			&u_l2g,
			&elementNeighborsArray,
			&u_dof,
			&u_dof_limit,
			&dm,
			&w,
			&v,
			&dV,
			&rhs,
			&elementResidual,
			&slumpedMassMatrixCorrection))
    return NULL;
  //1dv2 is mass conservative attempt
  calculateBerzinsSlumpedMassApproximation1dv2(SHAPE(u_l2g)[0],
					    SHAPE(elementNeighborsArray)[1],
					    SHAPE(dm)[1],
					    SHAPE(v)[2],
					    SHAPE(w)[2],
					    IDATA(u_l2g),
					    IDATA(elementNeighborsArray),
					    DDATA(u_dof),
					    DDATA(u_dof_limit),
					    DDATA(dm),
					    DDATA(w),
					    DDATA(v),
					    DDATA(dV),
					    DDATA(rhs),
					    DDATA(elementResidual),
					    DDATA(slumpedMassMatrixCorrection));

  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_calculateBerzinsSlumpedMassApproximation2d(PyObject* self,
								  PyObject* args)
{
  PyObject *u_l2g,
    *elementNeighborsArray,
    *u_dof,*u_dof_limit,
    *dm,*w,*v,
    *dV,*rhs,
    *slumpedMassMatrixCorrection,
    *elementResidual; 

  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOO",
			&u_l2g,
			&elementNeighborsArray,
			&u_dof,
			&u_dof_limit,
			&dm,
			&w,
			&v,
			&dV,
			&rhs,
			&elementResidual,
			&slumpedMassMatrixCorrection))
    return NULL;
  
  calculateBerzinsSlumpedMassApproximation2d(SHAPE(u_l2g)[0],
					    SHAPE(elementNeighborsArray)[1],
					    SHAPE(dm)[1],
					    SHAPE(v)[2],
					    SHAPE(w)[2],
					    IDATA(u_l2g),
					    IDATA(elementNeighborsArray),
					    DDATA(u_dof),
					    DDATA(u_dof_limit),
					    DDATA(dm),
					    DDATA(w),
					    DDATA(v),
					    DDATA(dV),
					    DDATA(rhs),
					    DDATA(elementResidual),
					    DDATA(slumpedMassMatrixCorrection));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cellam_updateElementJacobianWithSlumpedMassCorrection(PyObject* self,
								       PyObject* args)
{
  PyObject *elementMassMatrixCorrection,
    *elementJacobian; 

  if (!PyArg_ParseTuple(args,
                        "OO",
			&elementMassMatrixCorrection,
			&elementJacobian))
    return NULL;
  
  updateElementJacobianWithSlumpedMassCorrection(SHAPE(elementJacobian)[0],
						 SHAPE(elementJacobian)[2],
						 SHAPE(elementJacobian)[1],
						 DDATA(elementMassMatrixCorrection),
						 DDATA(elementJacobian));

  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cellam_manuallyUpdateGlobalMassMatrix(PyObject* self,
						       PyObject* args)
{
  PyObject *rowptr,*colind,*u_l2g,
    *u_dof,
    *dm,*w,*v,
    *dV,
    *globalMassMatrix;

  if (!PyArg_ParseTuple(args,
			"OOOOOOOOO",
			&rowptr,
			&colind,
			&u_l2g,
			&u_dof,
			&dm,
			&w,
			&v,
			&dV,
			&globalMassMatrix))
    return NULL;
  
  manuallyUpdateGlobalMassMatrix(SHAPE(u_l2g)[0],
				 SHAPE(dm)[1],
				 SHAPE(v)[2],
				 SHAPE(w)[2],
				 IDATA(rowptr),
				 IDATA(colind),
				 IDATA(u_l2g),
				 DDATA(u_dof),
				 DDATA(dm),
				 DDATA(w),
				 DDATA(v),
				 DDATA(dV),
				 DDATA(globalMassMatrix));
  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cellam_calculateElementSlumpedMassApproximationFromGlobalEdgeLimiter(PyObject* self,
										      PyObject* args)
{
  PyObject *rowptr,*colind,*u_l2g,
    *u_dof,
    *dm,*w,*v,
    *dV,
    *globalEdgeSlumpingParameter,
    *slumpedMassMatrixCorrection,
    *elementResidual; 

  if (!PyArg_ParseTuple(args,
			"OOOOOOOOOOO",
			&rowptr,
			&colind,
			&u_l2g,
			&u_dof,
			&dm,
			&w,
			&v,
			&dV,
			&globalEdgeSlumpingParameter,
			&elementResidual,
			&slumpedMassMatrixCorrection))
    return NULL;
  
  calculateElementSlumpedMassApproximationFromGlobalEdgeLimiter(SHAPE(u_l2g)[0],
								SHAPE(dm)[1],
								SHAPE(v)[2],
								SHAPE(w)[2],
								IDATA(rowptr),
								IDATA(colind),
								IDATA(u_l2g),
								DDATA(u_dof),
								DDATA(dm),
								DDATA(w),
								DDATA(v),
								DDATA(dV),
								DDATA(globalEdgeSlumpingParameter),
								DDATA(elementResidual),
								DDATA(slumpedMassMatrixCorrection));
  
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_computeSlumpingParametersFCT_KuzminMoeller10(PyObject* self,
								     PyObject* args)
{
  PyObject *rowptr,*colind,
    *u_dof,*u_dof_limit,
    *Mc,
    *Rip,*Rim,
    *globalEdgeSlumpingParameter;

  if (!PyArg_ParseTuple(args,
			"OOOOOOOO",
			&rowptr,
			&colind,
			&u_dof,
			&u_dof_limit,
			&Mc,
			&Rip,
			&Rim,
			&globalEdgeSlumpingParameter))
    return NULL;
  
  computeSlumpingParametersFCT_KuzminMoeller10(SHAPE(rowptr)[0]-1,
					       IDATA(rowptr),
					       IDATA(colind),
					       DDATA(u_dof),
					       DDATA(u_dof_limit),
					       DDATA(Mc),
					       DDATA(Rip),
					       DDATA(Rim),
					       DDATA(globalEdgeSlumpingParameter));
  
  Py_INCREF(Py_None); 
  return Py_None;
}

//-- SSIP stuff below --
static PyObject* cellam_generateQuadratureArraysForSSIPs(PyObject* self,
							 PyObject* args)

{
  //input 
  PyObject *nodeArray,
    *elementNodesArray,
    *elementBoundariesArray,
    *elementBoundaryLocalOuterNormalsArray,
    *elementBoundaryBarycentersArray,
    *element_track,
    *flag_track,
    *x_track,
    *q_x,
    *q_dV;
  double boundaryTolerance,
    neighborTolerance;
  //output
  PyObject *x_gq,*dV_gq,*elements_gq;
  int nPointsTracked = 1;
  if(!PyArg_ParseTuple(args,
                       "ddOOOOOOOOOO",
                       &boundaryTolerance,
		       &neighborTolerance,
		       &nodeArray,
		       &elementNodesArray,
		       &elementBoundariesArray,
		       &elementBoundaryLocalOuterNormalsArray,
		       &elementBoundaryBarycentersArray,
		       &element_track,
		       &flag_track,
		       &x_track,
		       &q_x,
		       &q_dV))
    return NULL;
  for (int i = 0; i < ND(element_track); i++)
    nPointsTracked *= SHAPE(element_track)[i];
  //since just working through algorithm,
  //use stl containers to build points and copy over
  //rather than try to reuse memory etc
  std::vector<int> elements_gq_tmp;
  std::vector<double> x_gq_tmp;
  std::vector<double> dV_gq_tmp; 

  generateQuadratureArraysForSSIPs(SHAPE(elementNodesArray)[0],
				   SHAPE(elementNodesArray)[1],
				   SHAPE(elementBoundariesArray)[1],
				   SHAPE(elementBoundaryLocalOuterNormalsArray)[2],
				   nPointsTracked,
				   SHAPE(q_dV)[1],
				   boundaryTolerance,
				   neighborTolerance,
				   DDATA(nodeArray),
				   IDATA(elementNodesArray),
				   IDATA(elementBoundariesArray),
				   DDATA(elementBoundaryLocalOuterNormalsArray),
				   DDATA(elementBoundaryBarycentersArray),
				   IDATA(element_track),
				   IDATA(flag_track),
				   DDATA(x_track),
				   DDATA(q_x),
				   DDATA(q_dV),
				   elements_gq_tmp,
				   x_gq_tmp,
				   dV_gq_tmp);

  npy_intp nPoints_global = elements_gq_tmp.size();
  npy_intp dims[2];

  dims[0] = nPoints_global; dims[1] = 3;
  elements_gq = PyArray_SimpleNew(1,dims,NPY_INT);
  dV_gq       = PyArray_SimpleNew(1,dims,NPY_DOUBLE);
  x_gq        = PyArray_SimpleNew(2,dims,NPY_DOUBLE);

  //manually copy
  int* elements_gq_p = IDATA(elements_gq);
  double* dV_gq_p    = DDATA(dV_gq);
  double* x_gq_p     = DDATA(x_gq);

  for (int k=0; k < nPoints_global; k++)
    {
      //std::cout<<"cellam generateArraysForSSIPs k= "<<k;
      elements_gq_p[k] = elements_gq_tmp[k];
      dV_gq_p[k]       = dV_gq_tmp[k];
      //std::cout<<" ele tmp= "<<elements_gq_tmp[k] <<" ele_p= "<<elements_gq_p[k];
      //std::cout<<" dV tmp= "<<dV_gq_tmp[k] <<" dV_p= "<<dV_gq_p[k];
      for (int I=0; I < 3; I++)
  	{
  	  x_gq_p[k*3+I]        = x_gq_tmp[k*3+I];
  	  //std::cout<<" I= "<<I<<"  x tmp= "<<x_gq_tmp[k*3+I] <<" x_p= "<<x_gq_p[k*3+I];
  	}
      //std::cout<<std::endl;
    }

  return Py_BuildValue("OOO",
		       elements_gq,
		       dV_gq,
		       x_gq);
}


static PyObject* cellam_generateArraysForTrackedSSIPs(PyObject* self,
						      PyObject* args)

{
  //input 
  PyObject *nodeArray,
    *elementNodesArray,
    *elementBoundariesArray,
    *elementBoundaryLocalOuterNormalsArray,
    *elementBoundaryBarycentersArray,
    *element_track,
    *flag_track,
    *x_track;
  double boundaryTolerance,
    neighborTolerance;
  //output
  PyObject *x_ssip,*element_offsets_ssip;
  int nPointsTracked = 1;
  if(!PyArg_ParseTuple(args,
                       "ddOOOOOOOO",
                       &boundaryTolerance,
		       &neighborTolerance,
		       &nodeArray,
		       &elementNodesArray,
		       &elementBoundariesArray,
		       &elementBoundaryLocalOuterNormalsArray,
		       &elementBoundaryBarycentersArray,
		       &element_track,
		       &flag_track,
		       &x_track))
    return NULL;
  for (int i = 0; i < ND(element_track); i++)
    nPointsTracked *= SHAPE(element_track)[i];
  //since just working through algorithm,
  //use stl containers to build points and copy over
  //rather than try to reuse memory etc
  
  std::vector<int> element_offsets_ssip_tmp(SHAPE(elementNodesArray)[0]+1);
  std::vector<double> x_ssip_tmp;

  generateArraysForTrackedSSIPs(SHAPE(elementNodesArray)[0],
				SHAPE(elementNodesArray)[1],
				SHAPE(elementBoundariesArray)[1],
				SHAPE(elementBoundaryLocalOuterNormalsArray)[2],
				nPointsTracked,
				boundaryTolerance,
				neighborTolerance,
				DDATA(nodeArray),
				IDATA(elementNodesArray),
				IDATA(elementBoundariesArray),
				DDATA(elementBoundaryLocalOuterNormalsArray),
				DDATA(elementBoundaryBarycentersArray),
				IDATA(element_track),
				IDATA(flag_track),
				DDATA(x_track),
				element_offsets_ssip_tmp,
				x_ssip_tmp);

  npy_intp nPoints_global = x_ssip_tmp.size()/3;
  npy_intp nElement_offsets = element_offsets_ssip_tmp.size();
  npy_intp dims[2];

  dims[0] = nElement_offsets; dims[1] = 0;
  element_offsets_ssip = PyArray_SimpleNew(1,dims,NPY_INT);
  dims[0] = nPoints_global; dims[1] = 3;
  x_ssip     = PyArray_SimpleNew(2,dims,NPY_DOUBLE);

  //manually copy
  int* element_offsets_ssip_p = IDATA(element_offsets_ssip);
  double* x_ssip_p     = DDATA(x_ssip);

  for (int k = 0; k < nElement_offsets; k++)
    {
      //std::cout<<"cellam generateArraysForSSIPs k= "<<k;
      element_offsets_ssip_p[k] = element_offsets_ssip_tmp[k];
    }
  for (int k=0; k < nPoints_global; k++)
    {
      //std::cout<<"cellam generateArraysForSSIPs k= "<<k;
      for (int I=0; I < 3; I++)
  	{
  	  x_ssip_p[k*3+I]        = x_ssip_tmp[k*3+I];
  	  //std::cout<<" I= "<<I<<"  x tmp= "<<x_gq_tmp[k*3+I] <<" x_p= "<<x_gq_p[k*3+I];
  	}
      //std::cout<<std::endl;
    }

  return Py_BuildValue("OO",
		       element_offsets_ssip,
		       x_ssip);
}
static PyObject* cellam_accumulateSourceContribution(PyObject* self,
						     PyObject* args)
{
  //input
  PyObject * traj_offsets, //traj_offsets[i] = start of trajectory info for particle i
		           //n_i = traj_offsets[i+1]-traj_offsets[i] 
     * x_traj, //particle trajectories: x (size 3\sum_i n_i])
     * t_traj, //particle trajectories: t
     * elem_traj, //particle trajectories: element id's
     * massSource_t, //discrete t values (knot's) for mass source
     * massSource_m, //values for mass source at knot's
     * decay_coef_element,  //linear decay: nParticleFlags * nElements_global
     * retardation_factor_element,  //Retardation factor: nParticleFlags * nElements_global
     * particleFlags, //The particle 'type' associated with particles
     * c_element; //element concentrations
  double tau;
  if(!PyArg_ParseTuple(args,
                       "dOOOOOOOOOO",
		       &tau,
                       &traj_offsets,
		       &x_traj,
		       &t_traj,
		       &elem_traj,
		       &massSource_t,
		       &massSource_m,
		       &decay_coef_element,
		       &retardation_factor_element,
		       &particleFlags,
		       &c_element))
    return NULL;
 
  accumulateSourceContribution(SHAPE(traj_offsets)[0]-1,
			       SHAPE(c_element)[0],
			       SHAPE(decay_coef_element)[1],
			       SHAPE(massSource_t)[0],
			       tau,           
			       IDATA(traj_offsets),
			       DDATA(x_traj),
			       DDATA(t_traj),
			       IDATA(elem_traj),
			       DDATA(massSource_t),
			       DDATA(massSource_m),
			       DDATA(decay_coef_element),
			       DDATA(retardation_factor_element),
			       IDATA(particleFlags),
			       DDATA(c_element));

  
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cellam_accumulateSourceContributionMaterialTypes(PyObject* self,
								  PyObject* args)
{
  //input
  PyObject * traj_offsets, //traj_offsets[i] = start of trajectory info for particle i
		           //n_i = traj_offsets[i+1]-traj_offsets[i] 
     * x_traj, //particle trajectories: x (size 3\sum_i n_i])
     * t_traj, //particle trajectories: t
     * elem_traj, //particle trajectories: element id's
     * massSource_t, //discrete t values (knot's) for mass source
     * massSource_m, //values for mass source at knot's
     * material_types_element, //identifier for material of each element 
     * decay_coef_types,  //linear decay: nMaterialTypes * nParticleFlags  
     * retardation_factor_types,  //Retardation factor: nMaterialTypes * nParticleFlags 
     * particle_flags, //The particle 'type' associated with particles
     * c_element; //element concentrations
  double tau;
  if(!PyArg_ParseTuple(args,
                       "dOOOOOOOOOOO",
		       &tau,
                       &traj_offsets,
		       &x_traj,
		       &t_traj,
		       &elem_traj,
		       &massSource_t,
		       &massSource_m,
		       &material_types_element,
		       &decay_coef_types,
		       &retardation_factor_types,
		       &particle_flags,
		       &c_element))
    return NULL;
 
  accumulateSourceContributionMaterialTypes(SHAPE(traj_offsets)[0]-1,
					    SHAPE(c_element)[0],
					    SHAPE(decay_coef_types)[1],
					    SHAPE(massSource_t)[0],
					    tau,           
					    IDATA(traj_offsets),
					    DDATA(x_traj),
					    DDATA(t_traj),
					    IDATA(elem_traj),
					    DDATA(massSource_t),
					    DDATA(massSource_m),
					    IDATA(material_types_element),
					    DDATA(decay_coef_types),
					    DDATA(retardation_factor_types),
					    IDATA(particle_flags),
					    DDATA(c_element));

  
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cellam_integratePiecewiseLinearMassSource(PyObject* self,
							   PyObject* args)
{
  //input
  PyObject * t_vals,
    * m_vals;
  double t_in,t_out,tau_out,tau_in,decay;

  if(!PyArg_ParseTuple(args,
                       "dddddOO",
		       &t_in,
		       &t_out,
		       &tau_out,
		       &tau_in,
		       &decay,
		       &t_vals,
		       &m_vals))
    return NULL;
 
  double massint = integratePiecewiseLinearMassSource(SHAPE(t_vals)[0],
						     DDATA(t_vals),
						     DDATA(m_vals),
						     t_in,
						     t_out,
						     tau_out,
						     tau_in,
						     decay);
  
  return Py_BuildValue("d",
		       massint);

}
static PyObject* cellam_volume123(PyObject* self,
				  PyObject* args)
{
  //input
  PyObject * elementNodes_element,
    * nodeArray;
  int nSpace;

  if(!PyArg_ParseTuple(args,
                       "iOO",
		       &nSpace,
		       &elementNodes_element,
		       &nodeArray))
    return NULL;
 
  double volume = volume123(nSpace,
			    SHAPE(elementNodes_element)[0],
			    IDATA(elementNodes_element),
			    DDATA(nodeArray));
  return Py_BuildValue("d",
		       volume);

}
static PyObject* cellam_nodalProjection123(PyObject* self,
					   PyObject* args)
{
  //input
  PyObject * elementNodes_offset,
    * elementNodesArray,
    * elementVolumesArray,
    * c_ele,
    * nodeStar_volume,
    * c_nodal;
  int calculateNodeStarVolume;
  if(!PyArg_ParseTuple(args,
                       "iOOOOOO",
		       &calculateNodeStarVolume,
		       &elementNodes_offset,
		       &elementNodesArray,
		       &elementVolumesArray,
		       &c_ele,
		       &nodeStar_volume,
		       &c_nodal))
    return NULL;
  nodalProjection123(calculateNodeStarVolume,
		     SHAPE(c_ele)[0],
		     SHAPE(c_nodal)[0],
		     IDATA(elementNodes_offset),
		     IDATA(elementNodesArray),
		     DDATA(elementVolumesArray),
		     DDATA(c_ele),
		     DDATA(nodeStar_volume),
		     DDATA(c_nodal));

  Py_INCREF(Py_None); 
  return Py_None;

}
static PyMethodDef cellamMethods[] = {
 { "updateOldMass_weak",
    cellam_updateOldMass_weak,
   METH_VARARGS, 
   "update ellam element residual with mass accumulation from previous time level"},
 { "updateOldMass_weak_arbitraryQuadrature",
    cellam_updateOldMass_weak_arbitraryQuadrature,
   METH_VARARGS, 
   "update ellam element residual with mass accumulation from previous time level using genera set of quadrature points"},
 { "updateNewMass_weak",
    cellam_updateNewMass_weak,
   METH_VARARGS, 
   "update ellam element residual with mass accumulation from current time level"},
 { "evaluateSolutionAtTrackedPoints",
    cellam_evaluateSolutionAtTrackedPoints,
   METH_VARARGS, 
   "evaluate solution at tracked points"},
 { "updateExteriorOutflowBoundaryFlux",
    cellam_updateExteriorOutflowBoundaryFlux,
   METH_VARARGS, 
   "approximate outflow boundary integral using trapezoidal rule"},
 { "updateExteriorOutflowBoundaryFluxInGlobalResidual",
    cellam_updateExteriorOutflowBoundaryFluxInGlobalResidual,
   METH_VARARGS, 
   "approximate outflow boundary integral using trapezoidal rule directly in global residual"},
 { "updateExteriorOutflowBoundaryFluxJacobian",
   cellam_updateExteriorOutflowBoundaryFluxJacobian,
   METH_VARARGS, 
   "update jacobian to reflext approximate outflow boundary integral using trapezoidal rule"},
 { "updateExteriorOutflowBoundaryFluxGlobalJacobian",
   cellam_updateExteriorOutflowBoundaryFluxGlobalJacobian,
   METH_VARARGS, 
   "update global jacobian to reflext approximate outflow boundary integral using trapezoidal rule"},
 { "markInflowBoundaryPoints",
   cellam_markInflowBoundaryPoints,
   METH_VARARGS, 
   "mark points that need to be integrated for inflow boundary"},
 { "accumulateInflowFlux",
   cellam_accumulateInflowFlux,
   METH_VARARGS, 
   "apply inflow boundary contribution to local element residual "},
 { "accumulateInflowFluxInGlobalResidual",
   cellam_accumulateInflowFluxInGlobalResidual,
   METH_VARARGS, 
   "apply inflow boundary contribution to residual "},
 { "tagNegligibleIntegrationPoints",
   cellam_tagNegligibleIntegrationPoints,
   METH_VARARGS, 
   "tag integration points with magnitude less than tolerance "},
 { "calculateSlumpedMassApproximation1d",
    cellam_calculateSlumpedMassApproximation1d,
   METH_VARARGS, 
   "calculate 1d slumping approximation from Russell and Binning and apply to residual"},
 { "calculateSlumpedMassApproximation1d_local",
    cellam_calculateSlumpedMassApproximation1d_local,
   METH_VARARGS, 
   "calculate 1d slumping approximation from Russell and Binning just looking at local element condition and apply to residual"},
 { "calculateSlumpedMassApproximation2d",
    cellam_calculateSlumpedMassApproximation2d,
   METH_VARARGS, 
   "calculate 2d slumping approximation for just local mass matrix and apply to residual"},
 { "calculateSlumpedMassApproximation2d_upwind",
    cellam_calculateSlumpedMassApproximation2d_upwind,
   METH_VARARGS, 
   "calculate 2d slumping approximation following 1d paradigm in upwind direction and apply to residual"},
 { "calculateBerzinsSlumpedMassApproximation1d",
    cellam_calculateBerzinsSlumpedMassApproximation1d,
   METH_VARARGS, 
   "calculate 1d slumping approximation from Berzins and apply to residual"},
 { "calculateBerzinsSlumpedMassApproximation2d",
    cellam_calculateBerzinsSlumpedMassApproximation2d,
   METH_VARARGS, 
   "calculate 2d slumping approximation from Berzins and apply to residual"},
 { "updateElementJacobianWithSlumpedMassApproximation",
    cellam_updateElementJacobianWithSlumpedMassApproximation,
   METH_VARARGS, 
   "update jacobian matrix to account for slumping"},
 { "updateElementJacobianWithSlumpedMassCorrection",
    cellam_updateElementJacobianWithSlumpedMassCorrection,
   METH_VARARGS, 
   "update jacobian matrix to account for slumping"},
 { "manuallyUpdateGlobalMassMatrix",
   cellam_manuallyUpdateGlobalMassMatrix,
   METH_VARARGS, 
   "update a global mass matrix approx for all dofs "},
 { "calculateElementSlumpedMassApproximationFromGlobalEdgeLimiter",
   cellam_calculateElementSlumpedMassApproximationFromGlobalEdgeLimiter,
   METH_VARARGS, 
   "update element residual and jacobian based on FCT limiting parameters computed globally to keep assembly structure for now "},
 { "computeSlumpingParametersFCT_KuzminMoeller10",
   cellam_computeSlumpingParametersFCT_KuzminMoeller10,
   METH_VARARGS, 
   "compute slumping parameters based on FCT approach from Kuzmin Moeller10 assuming have an assembled global mass matrix"},
 //SSIPs
 { "generateQuadratureArraysForSSIPs",
   cellam_generateQuadratureArraysForSSIPs,
   METH_VARARGS, 
   "return point to element look up array, quadrature weight array, and quadrature point array for SSIPs"},
 { "generateArraysForTrackedSSIPs",
   cellam_generateArraysForTrackedSSIPs,
   METH_VARARGS, 
   "return element offsets array and point array for subset of SSIPs that actually get tracked to an element"},
 { "accumulateSourceContribution",
   cellam_accumulateSourceContribution,
   METH_VARARGS, 
   "apply convolution-based particle tracking (CBPT) approach to generate element-based concentrations"},
 { "accumulateSourceContributionMaterialTypes",
   cellam_accumulateSourceContributionMaterialTypes,
   METH_VARARGS, 
   "apply convolution-based particle tracking (CBPT) approach to generate element-based concentrations assuming material types assigned using element material id and particle id"},
 { "integratePiecewiseLinearMassSource",
   cellam_integratePiecewiseLinearMassSource,
   METH_VARARGS, 
   "integrate piecewise-linear mass source with exponential decay for CBPT"},
 { "volume123",
   cellam_volume123,
   METH_VARARGS, 
   "compute volume for simplex, quadrilateral, triangular prism, or hexahedron"},
 { "nodalProjection123",
   cellam_nodalProjection123,
   METH_VARARGS, 
   "mass-lumped L2 projection for simplex, quadrilateral, triangular prism, or hexahedron"},
 { NULL,NULL,0,NULL}
};

extern "C"
{
PyMODINIT_FUNC initcellam(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cellam", cellamMethods);
  d = PyModule_GetDict(m);
  import_array();
}
}//extern "C"
