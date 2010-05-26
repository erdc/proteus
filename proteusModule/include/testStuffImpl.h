#ifndef TESTSTUFFIMPL_H
#define TESTSTUFFIMPL_H
#include "Python.h"
#include PROTEUS_LAPACK_H

/** 
 \file testStuffImpl.h
 \brief Pyton interface to things that don't yet  have a home
*/

/**
   \defgroup testStuffImpl testStuffImpl
   \brief Python interface to things that don't yet have a home
   @{ 
*/


/***********************************************************************
  Apply Algorithm 5.2 in Barth and Sethian 98 for level set equation

  nElements_global : IN, number of elements (Ne)
  nDOF_element     : IN, number of dofs per element := nSpace+1
  nSpace           : IN, space dimension (Nsd)
  nQuadraturePoints_element: IN, number of quadrature points per element (nq)
  nQuadraturePoints_elementBoundary : IN, number of quadrature points per
                    element boundary (nqb)

  elementDiameters : IN, Ne x 1. 
                     elementDiameters[I]=h_{E_I} for global element I
  elementUnitOuterNormals : IN, Ne x Nsd+1 x nqb x Nsd. unit outer normals for element boundaries
                     elementUnitOuterNormals[I,j,k,:] = \vec n^j_I(x^j_k)
  detJ             : IN, Ne x nq.  
                     detJ[I,k] = |det J_{I,k}| for global element I and 
                     quadrature point k
                     = |E_I|/[1,2,6] for 1d,2d,3d
  sqrtDetG         : IN, Ne x Nsd+1 x nqb. sqrt of metric tensor deterimant
                     for local boundary j of global element I eval'd at 
                     quadrature point k
                     sqrtDetG[I,j,k] = sqrt(det(g_I^j)(x_k))  
                     = |e_I^j|/[1,1,2]  for 1d,2d,3d
  elementQuadratureWeights : IN, nq x 1. 
                     element quadrature points on reference element

  l2g              : IN, Ne x ndof. local to global degree of freedom map
                     l2g[I,j] = global dof for element I's j'th degree of freedom
  phiIn            : IN, Ndof x 1. input degrees of freeom.
 
  H                : IN, Ne x nq. Hamiltonian F|\grad phi| eval'd at element
                     quad points
                     H[I,k] = H(\vec x^k_I)

  dH               : IN, Ne x nq x Nsd. derivative of H w.r.t \grad phi.
                     dH[I,k,:]=\pd{H}{\grad \phi}(\vec x^k_{I})

  r                : IN, Ne x nq. minus source term. 
                     r[I,k]  = -f(x^k_{I})
  dt               : IN, \Delta t
  phiOut           : IN, Ndof x 1. output degrees of freeom.
  wOut             : IN, Ndof x 1. output test function weights
  
 ***********************************************************************/
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
} NodeStarFactor;

extern void advanceStageP1_C0_GLS_lump(
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  int nQuadraturePoints_elementBoundary,
  double *elementDiameters,
  double *elementUnitOuterNormals,
  double *detJ,
  double *sqrtDetG,
  double *elementQuadratureWeights,
  int *l2g,
  double *phiIn,
  double *H,
  double *dH,
  double *r,
  double dt,
  double *phiOut,
  double *wOut
);
extern void advanceStageP1_C0_GLS_lump_noSource(
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  int nQuadraturePoints_elementBoundary,
  double *elementDiameters,
  double *elementUnitOuterNormals,
  double *detJ,
  double *sqrtDetG,
  double *elementQuadratureWeights,
  int *l2g,
  double *phiIn,
  double *H,
  double *dH,
  double dt,
  double *phiOut,
  double *wOut
);
extern void advanceStageP1_C0_SGS_lump_noSource(
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  int nQuadraturePoints_elementBoundary,
  double *elementDiameters,
  double *elementUnitOuterNormals,
  double *detJ,
  double *sqrtDetG,
  double *elementQuadratureWeights,
  int *l2g,
  double *phiIn,
  double *H,
  double *dH,
  double dt,
  double *phiOut,
  double *wOut
);
extern void advanceStageP1_C0_SUPG_lump(
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  int nQuadraturePoints_elementBoundary,
  double *elementDiameters,
  double *elementUnitOuterNormals,
  double *detJ,
  double *sqrtDetG,
  double *elementQuadratureWeights,
  int *l2g,
  double *phiIn,
  double *H,
  double *dH,
  double *r,
  double dt,
  double *phiOut,
  double *wOut
);
extern void advanceStageRedistanceP1_C0_SUPG_lump(
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  int nQuadraturePoints_elementBoundary,
  double *elementDiameters,
  double *elementUnitOuterNormals,
  double *detJ,
  double *sqrtDetG,
  double *elementQuadratureWeights,
  int *l2g,
  double *phiIn,
  double *phiEvalIn,
  double dt,
  double eps,
  double *phiOut,
  double *wOut
);
extern void advanceStageRedistanceWeakDirP1_C0_SUPG_lump(
  int nElements_global,
  int nElementBoundaries_element,
  int nDOF_element,
  int nSpace,
  int nQuadraturePoints_element,
  int nQuadraturePoints_elementBoundary,
  int nExteriorElementBoundaries_global,
  double *elementDiameters,
  double *elementUnitOuterNormals,
  double *detJ,
  double *sqrtDetG,
  double *elementQuadratureWeights,
  int *l2g,
  int *exteriorElementBoundariesArray,
  int *elementBoundaryElementsArray,
  double *phiIn,
  double *phiEvalIn,
  double dt,
  double eps,
  int *weakDirichletFlag,
  double *phiOut,
  double *wOut
);
extern void calculateConservationResidualPWLv2(
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
  int *nodeStarOffsets,
  int *nElements_node,
  int *fluxElementBoundaries,
  double *elementResidual,
  double *vAverage,
  double *starU,
  double *dX,
  double *w,
  double *normal,
  double *conservationResidual,
  double *starR,
  double *vConservative,
  double *vConservative_element
);
extern void calculateConservationJacobianPWLv2(
  int nNodes_global,
  int nNodes_internal,
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nNodes_element,
  int nSpace,
  int nComponentsTotal,
  int componentId,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElements,
  int *elementBoundaryLocalElementBoundaries,
  int *elementNodes,
  int *nodeStarElements,
  int *nodeStarElementNeighbors,
  int *nodeStarOffsets,
  int *nodeStarJacobianOffsets,
  int *nElements_node,
  int *internalNodes,
  int *fluxElementBoundaries,
  int *fluxBoundaryNodes,
  double *w,
  double *normal,
  double *starJacobian
);
extern void calculateConservationFluxPWLv2(
  int nNodes_global,
  int nNodes_internal,
  int componentId,
  int *nElements_node,
  int *nodeStarOffsets,
  int *nodeStarJacobianOffsets,
  int *internalNodes,
  int *fluxBoundaryNodes,
  double *starR,
  double *starJ,
  double *starU
);
extern void loadBoundaryFluxIntoGlobalElementBoundaryVelocity(
  int nExteriorElementBoundaries_global,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *exteriorElementBoundaries,
  int *fluxElementBoundaries,
  double *normal,
  double *flux,
  double updateCoef,
  double *velocity
);
extern void copyGlobalElementBoundaryVelocityToElementBoundary(
  int nElements_global,
  int nInteriorElementBoundaries_global,
  int nExteriorElementBoundaries_global,
  int nElementBoundaries_global,
  int nElementBoundaries_element,
  int nQuadraturePoints_elementBoundary,
  int nSpace,
  int *interiorElementBoundaries,
  int *exteriorElementBoundaries,
  int *elementBoundaryElementsArray,
  int *elementBoundaryLocalElementBoundariesArray,
  double *velocityBoundary_global,
  double *velocityBoundary_element
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
  PROTEUS_LAPACK_INTEGER *** subdomain_pivots_p
);
extern int nodeStar_free(
  int N,
  int *subdomain_dim,
  double **subdomain_L,
  double **subdomain_R,
  double **subdomain_U,
  PROTEUS_LAPACK_INTEGER ** subdomain_pivots
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
  int *N_p,
  int **subdomain_dim_p,
  double ***subdomain_L_p,
  double ***subdomain_R_p,
  double ***subdomain_U_p,
  PROTEUS_LAPACK_INTEGER *** subdomain_pivots_p
);
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
extern void calculateConservationFluxPWLv3(
  int nNodes_global,
  int nNodes_internal,
  int *nElements_node,
  int *internalNodes,
  int *fluxBoundaryNodes,
  NodeStarFactor * nodeStarFactor
);
/** @} */
#endif
