#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H
#include "Python.h"
#include PROTEUS_LAPACK_H

/*!
 \file postprocessing.h
 \brief Python interface to velocity postprocessing library.
*/

/**
   \defgroup postprocessing postprocessing
   \brief Python interface to velocity postprocessing library.
   @{ 
*/

/***********************************************************************
  data structure for node star solves, same as ASM smoother basically
 ***********************************************************************/
typedef struct 
{
  PyObject_HEAD
  int N;
  int *subdomain_dim;
  /*int **l2g_L;*/
  double **subdomain_L,
    **subdomain_R,
    **subdomain_U;
  PROTEUS_LAPACK_INTEGER** subdomain_pivots;
  PROTEUS_LAPACK_INTEGER** subdomain_column_pivots;
} NodeStarFactor;

extern void invertLocal(
  int nSpace,
  double (*A)[3],
  double (*AI)[3]
);
extern void updateSelectedExteriorElementBoundaryFlux(int nExteriorElementBoundaries_global,
					       int nElementBoundaries_element,
					       int nQuadraturePoints_elementBoundary,
					       int nDOF_test_element,
					       int* exteriorElementBoundaries,
					       int* elementBoundaryElements,
					       int* elementBoundaryLocalElementBoundaries,
					       int* skipflag_elementBoundaries,
					       double* flux,
					       double* w_dS,
					       double* residual);
extern void updateRT0velocityWithAveragedPotentialP1nc(int nElements_global,
						int nQuadraturePoints_element,
						int nSpace,
						double * detJ,
						double * quad_a,
						double * phi,
						double * gradphi,
						double * a, 
						double * rt0vdofs);
extern void postProcessRT0velocityFromP1nc(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int * nFreeDOF_element,
  int * freeLocal_element,
  double * detJ,
  double * sqrt_det_g,
  double * n,
  double * elementBarycenters,
  double * quad_a,
  double * quad_f,
  double * w_dV_r,
  double * w_dV_m,
  double * u,
  double * gradu,
  double * a, 
  double * f, 
  double * r, 
  double * mt,
  double * rt0vdofs
);
extern void postProcessRT0velocityFromP1ncNoMass(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int * nFreeDOF_element,
  int * freeLocal_element,
  double * detJ,
  double * sqrt_det_g,
  double * n,
  double * elementBarycenters,
  double * quad_a,
  double * quad_f,
  double * w_dV_r,
  double * u,
  double * gradu,
  double * a, 
  double * f, 
  double * r, 
  double * rt0vdofs
);
extern void postProcessRT0velocityFromP1ncV2(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *nFreeDOF_element,
  int *freeLocal_element,
  double *detJ,
  double *sqrt_det_g,
  double *n,
  double *elementBarycenters,
  double *quad_a,
  double *quad_f,
  double *w_dV_r,
  double *w_dV_m,
  double *u,
  double *gradu,
  double *a,
  double *f,
  double *r,
  double *mt,
  double *rt0vdofs
);
extern void postProcessRT0velocityFromP1ncV2noMass(
  int nElements_global,
  int nQuadraturePoints_element,
  int nDOF_test_element,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *nFreeDOF_element,
  int *freeLocal_element,
  double *detJ,
  double *sqrt_det_g,
  double *n,
  double *elementBarycenters,
  double *quad_a,
  double *quad_f,
  double *w_dV_r,
  double *u,
  double *gradu,
  double *a,
  double *f,
  double *r,
  double *rt0vdofs
);
extern void getElementRT0velocityValues(
  int nElements_global,
  int nPoints_element,
  int nSpace,
  double *x_element,
  double *rt0vdofs_element,
  double *v_element
);
extern void getElementBoundaryRT0velocityValues(
  int nElements_global,
  int nElementBoundaries_element,
  int nPoints_elementBoundary,
  int nSpace,
  double *x_elementBoundary,
  double *rt0vdofs_element,
  double *v_elementBoundary
);
extern void getGlobalElementBoundaryRT0velocityValues(
  int nElementBoundaries_global,
  int nPoints_elementBoundary,
  int nSpace,
  int *elementBoundaryElementsArray,
  double *x_elementBoundary_global,
  double *rt0vdofs_element,
  double *v_elementBoundary_global
);
extern void postProcessRT0potentialFromP1nc(
  int nElements_global,
  int nQuadraturePoints_element,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  double *uQuadratureWeights_element,
  double *elementBarycenters,
  double *aElementQuadratureWeights,
  double *detJ,
  double *uQuadratureWeights_elementBoundary,
  double *x,
  double *u,
  double *gradu,
  double *x_elementBoundary,
  double *u_elementBoundary,
  double *n,
  double *a,
  double *f,
  double *r,
  double *rt0vdofs,
  double *rt0potential
);
extern void projectElementBoundaryVelocityToRT0fluxRep(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  double *elementBoundaryQuadratureWeights,
  double *n,
  double *v_elementBoundary,
  double *rt0vdofs_element
);
extern void projectElementBoundaryFluxToRT0fluxRep(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOF_RT0V_element,
  int *elementBoundaryElementsArray,
  int *elementBoundariesArray,
  double *elementBoundaryQuadratureWeights,
  double *flux_elementBoundary,
  double *rt0vdofs_element
);
extern void getElementRT0velocityValuesFluxRep(
  int nElements_global,
  int nElementBoundaries_element,
  int nPoints_element,
  int nSpace,
  int nDetVals_element,
  double *nodeArray,
  int *elementNodesArray,
  double *abs_det_J,
  double *x_element,
  double *rt0vdofs_element,
  double *v_element
);
extern void getElementBoundaryRT0velocityValuesFluxRep(
  int nElements_global,
  int nElementBoundaries_element,
  int nPoints_elementBoundary,
  int nSpace,
  int nDetVals_element,
  double *nodeArray,
  int *elementNodesArray,
  double *abs_det_J,
  double *x_elementBoundary,
  double *rt0vdofs_element,
  double *v_elementBoundary
);
extern void getGlobalElementBoundaryRT0velocityValuesFluxRep(
  int nElementBoundaries_global,
  int nPoints_elementBoundary_global,
  int nSpace,
  int nDetVals_element,
  double *nodeArray,
  int *elementNodesArray,
  int *elementBoundaryElementsArray,
  double *abs_det_J,
  double *x_elementBoundary_global,
  double *rt0vdofs_element,
  double *v_elementBoundary_global
);
extern void buildLocalBDM1projectionMatrices(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int nDOFs_test_element,
  int nDOFs_trial_element,
  int nVDOFs_element,
  double *w_dS_f,
  double *ebq_n,
  double *ebq_v,
  double *BDMprojectionMat_element
);
extern void factorLocalBDM1projectionMatrices(
  int nElements_global,
  int nVDOFs_element,
  double *BDMprojectionMat_element,
  int *BDMprojectionMatPivots_element
);
extern void solveLocalBDM1projection(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int nDOFs_test_element,
  int nVDOFs_element,
  double *BDMprojectionMatFact_element,
  int *BDMprojectionMatPivots_element,
  double *w_dS_f,
  double *ebq_n,
  double *ebq_velocity,
  double *p1_velocity_dofs
);
extern void solveLocalBDM1projectionFromFlux(
  int nElements_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nDOFs_test_element,
  int nVDOFs_element,
  double *BDMprojectionMatFact_element,
  int *BDMprojectionMatPivots_element,
  int *elementBoundaryElementsArray,
  int *elementBoundariesArray,
  double *w_dS_f,
  double *ebq_global_flux,
  double *p1_velocity_dofs
);
extern void getElementBDM1velocityValuesLagrangeRep(
  int nElements_global,
  int nQuadraturePoints_element,
  int nSpace,
  int nDOF_trial_element,
  int nVDOF_element,
  double *q_v,
  double *p1_velocity_dofs,
  double *q_velocity
);
extern void calculateConservationResidualGlobalBoundaries(
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nNodes_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *exteriorElementBoundariesToSkip,
  double *dS,
  double *normal,
  double *elementResidual,
  double *velocity,
  double *conservationResidual
);
extern void sunWheelerGSsweep(
  int nElements_global,
  int nElementBoundaries_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *dS,
  double *normal,
  double *sqrt_det_g,
  double *alpha,
  double *fluxCorrection,
  double *conservationResidual
);
extern void fluxCorrectionVelocityUpdate(
  int nElements_global,
  int nElementBoundaries_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  double *dS,
  double *normal,
  double *fluxCorrection,
  double *vConservative,
  double *vConservative_element
);
extern void computeFluxCorrectionPWC(
  int nElementBoundaries_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  double *pwcW,
  double *pwcV,
  double *fluxCorrection
);
extern int nodeStar_init(
  int nElements_global,
  int nNodes_element,
  int nNodes_global,
  int *nElements_node,
  int *nodeStarElementsArray,
  int *nodeStarElementNeighborsArray,
  int *N_p,
  int **subdomain_dim_p,
  double ***subdomain_L_p,
  double ***subdomain_R_p,
  double ***subdomain_U_p,
  PROTEUS_LAPACK_INTEGER *** subdomain_pivots_p,
  PROTEUS_LAPACK_INTEGER *** subdomain_column_pivots_p
);
extern int nodeStar_free(
  int N,
  int *subdomain_dim,
  double **subdomain_L,
  double **subdomain_R,
  double **subdomain_U,
  PROTEUS_LAPACK_INTEGER ** subdomain_pivots,
  PROTEUS_LAPACK_INTEGER ** subdomain_column_pivots_p
);
extern int nodeStar_setU(
  NodeStarFactor * nodeStarFactor,
  double val
);
extern int nodeStar_copy(
  int other_N,
  int *other_subdomain_dim,
  double **other_subdomain_L,
  double **other_subdomain_R,
  double **other_subdomain_U,
  PROTEUS_LAPACK_INTEGER ** other_subdomain_pivots,
  PROTEUS_LAPACK_INTEGER ** other_subdomain_column_pivots,
  int *N_p,
  int **subdomain_dim_p,
  double ***subdomain_L_p,
  double ***subdomain_R_p,
  double ***subdomain_U_p,
  PROTEUS_LAPACK_INTEGER *** subdomain_pivots_p,
  PROTEUS_LAPACK_INTEGER *** subdomain_column_pivots_p
);
extern
void calculateConservationResidualPWL(int nElements_global,
				      int nInteriorElementBoundaries_global,
				      int nExteriorElementBoundaries_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_elementBoundary,
				      int nDOF_element,
				      int nNodes_element,
				      int nSpace,
				      int* interiorElementBoundaries,
				      int* exteriorElementBoundaries,
				      int* elementBoundaryElements,
				      int* elementBoundaryLocalElementBoundaries,
				      int* elementNodes,
				      int* dofMapl2g,
				      int* nodeStarElements,
				      int* nodeStarElementNeighbors,
				      int* nElements_node,
				      int* fluxElementBoundaries,
				      double* elementResidual,
				      double* vAverage,
				      double* dX,
				      double* w,
				      double* normal,
				      NodeStarFactor* nodeStarFactor,
				      double* conservationResidual,
				      double* vConservative,
				      double* vConservative_element);


extern
void calculateConservationJacobianPWL(int nNodes_global,
				      int nNodes_internal,
				      int nElements_global,
				      int nInteriorElementBoundaries_global,
				      int nExteriorElementBoundaries_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_elementBoundary,
				      int nNodes_element,
				      int nDOF_element,
				      int nSpace,
				      int* interiorElementBoundaries,
				      int* exteriorElementBoundaries,
				      int* elementBoundaryElements,
				      int* elementBoundaryLocalElementBoundaries,
				      int* elementNodes,
				      int* dofMapl2g,
				      int* nodeStarElements,
				      int* nodeStarElementNeighbors,
				      int* nElements_node,
				      int* internalNodes,
				      int* fluxElementBoundaries,
				      int* fluxBoundaryNodes,
				      double* w,
				      double* normal,
				      NodeStarFactor* nodeStarFactor);

extern
void calculateConservationFluxPWL(int nNodes_global,
				  int nNodes_internal,
				  int* nElements_node,
				  int* internalNodes,
				  int* fluxBoundaryNodes,
				  NodeStarFactor* nodeStarFactor);
extern
void calculateConservationResidualPWL_opt(int nNodes_owned,
				      int nElements_global,
				      int nInteriorElementBoundaries_global,
				      int nExteriorElementBoundaries_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_elementBoundary,
				      int nNodes_element,
				      int nSpace,
				      int* interiorElementBoundaries,
				      int* exteriorElementBoundaries,
				      int* elementBoundaryElements,
				      int* elementBoundaryLocalElementBoundaries,
				      int* elementNodes,
				      int* nodeStarElements,
				      int* nodeStarElementNeighbors,
				      int* nElements_node,
				      int* fluxElementBoundaries,
				      double* elementResidual,
				      double* vAverage,
				      double* dX,
				      double* w,
				      double* normal,
				      NodeStarFactor* nodeStarFactor,
				      double* conservationResidual,
				      double* vConservative,
				      double* vConservative_element);


extern
void calculateConservationJacobianPWL_opt(int nNodes_owned,
				      int nNodes_global,
				      int nNodes_internal,
				      int nElements_global,
				      int nInteriorElementBoundaries_global,
				      int nExteriorElementBoundaries_global,
				      int nElementBoundaries_element,
				      int nQuadraturePoints_elementBoundary,
				      int nNodes_element,
				      int nSpace,
				      int* interiorElementBoundaries,
				      int* exteriorElementBoundaries,
				      int* elementBoundaryElements,
				      int* elementBoundaryLocalElementBoundaries,
				      int* elementNodes,
				      int* nodeStarElements,
				      int* nodeStarElementNeighbors,
				      int* nElements_node,
				      int* internalNodes,
				      int* fluxElementBoundaries,
				      int* fluxBoundaryNodes,
				      double* w,
				      double* normal,
				      NodeStarFactor* nodeStarFactor);

extern
void calculateConservationFluxPWL_opt(int nNodes_owned,
				  int nNodes_global,
				  int nNodes_internal,
				  int* nElements_node,
				  int* internalNodes,
				  int* fluxBoundaryNodes,
				  NodeStarFactor* nodeStarFactor);

extern void calculateConservationResidualPWLv3(
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nNodes_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *elementNodes,
  int *nodeStarElements,
  int *nodeStarElementNeighbors,
  int *nElements_node,
  int *fluxElementBoundaries,
  double *elementResidual,
  double *vAverage,
  double *dX,
  double *w,
  double *normal,
  NodeStarFactor * nodeStarFactor,
  double *conservationResidual,
  double *vConservative,
  double *vConservative_element
);
extern void calculateConservationJacobianPWLv3(
  int nNodes_global,
  int nNodes_internal,
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nNodes_element,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *elementNodes,
  int *nodeStarElements,
  int *nodeStarElementNeighbors,
  int *nElements_node,
  int *internalNodes,
  int *fluxElementBoundaries,
  int *fluxBoundaryNodes,
  double *w,
  double *normal,
  NodeStarFactor * nodeStarFactor
);
extern void calculateConservationFluxPWLv3(int nNodes_global,
                                           int nNodes_internal,
                                           int* nElements_node,
                                           int* internalNodes,
                                           int* fluxBoundaryNodes,
                                           NodeStarFactor* nodeStarFactor);

extern void postprocessAdvectiveVelocityPointEval(int nPoints,
						  int nSpace,
						  double updateCoef,
						  const double* f,
						  double * velocity);

extern void postprocessDiffusiveVelocityPointEval(int nPoints,
						  int nSpace,
						  double updateCoef,
						  const double* a,
						  const double* grad_phi,
						  double * velocity);  
extern void calculateConservationResidualPWL_primative(int nElements_global,
						int nInteriorElementBoundaries_global,
						int nExteriorElementBoundaries_global,
						int nElementBoundaries_element,
						int nQuadraturePoints_elementBoundary,
						int nNodes_element,
						int nSpace,
						int* interiorElementBoundaries,
						int* exteriorElementBoundaries,
						int* elementBoundaryElements,
						int* elementBoundaryLocalElementBoundaries,
						int* skipflag_elementBoundaries,
						double* elementResidual,
						double* dX,
						double* normal,
						double* conservationResidual,
						double* vConservative);
/** @}*/
#endif
