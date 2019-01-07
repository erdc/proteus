#ifndef VTIMEINTEGRATION_H
#define VTIMEINTEGRATION_H
/*!
 \file timeIntegration.h
 \brief c interface for implementation of time discretization algorithms
*/

/*Pseudo-Transient Continuation TTE formula choice for dt*/
extern
void psiTCtteDT(int nPoints,
		double tau,
		double dtn,
		double dtnm1,
		double * yn,
		double * ypn,
		double * ypnm1,
		double * dtnp1);


/*************************************************************************
   DG limiting schemes
 *************************************************************************/
extern
void applyDGlimitingP1Lagrange1d(int limiterFlag,
				 int nElements_global,
				 int nNodes_element,
				 int nElementBoundaries_element,
				 int nDOF_element,
				 int * elementNodesArray,
				 int * elementNeighborsArray,
				 double * nodeArray,
				 double * elementBarycentersArray,
				 int * l2g,
				 int * tag,
				 double * Uin,
				 double * Uout);

extern
void applyDGlimitingP1Lagrange1d_withVacuumTol(int enforcePositivity,
	       double vacuumTol,
	       int nElements_global,
	       int nNodes_element,
	       int nElementBoundaries_element,
	       int nDOF_element,
	       int * elementNodesArray,
	       int * elementNeighborsArray,
	       double * nodeArray,
	       double * elementBarycentersArray,
	       int * l2g,
	       int * tag,
	       double * Uin,
	       double * Uout);
extern
void computeElementNeighborShapeGradients(int nElements_global,
					  int nElementBoundaries_element,
					  int nSpace,
					  const int * elementBoundariesArray,
					  const int * elementNeighborsArray,
					  double * elementBarycentersArray,
					  double * elementBoundaryBarycentersArray,
					  double * elementNeighborShapeGradients);
extern
void computeCockburnDGlimiterArrays2d(int nElements_global,
				      int nElementBoundaries_element,
				      int nSpace,
				      const int * elementBoundariesArray,
				      const int * elementNeighborsArray,
				      const double * elementBarycentersArray,
				      const double * elementBoundaryBarycentersArray,
				      const double * elementNeighborShapeGradients,
				      double * alphas,
				      int * alphaNeighbors);
extern
void applyCockburnDGlimiterP1Lagrange2d(double nu,
					double Mh2,
					int nElements_global,
					int nElementBoundaries_element,
					int nSpace,
					int nDOF_element,
					int * elementNeighborsArray,
					int * l2g,
					int * tag,
					double * alphas,
					int * alphaNeighbors,
					double * Uin,
					double * Uout);
extern
void applyDurlofskyDGlimiterP1Lagrange2d(int killExtrema,
					 int allowMinWithUndershoot,
					 int nElements_global,
					 int nElementBoundaries_element,
					 int nNodes_element,
					 int nSpace,
					 int nDOF_element,
					 const int * elementNeighborsArray,
					 const int * elementBoundariesArray,
					 const int * elementNodesArray,
					 const double * nodeArray,
					 const double * elementBarycentersArray,
					 const double * elementBoundaryBarycentersArray,
					 const double * elementNeighborShapeGradients,
					 const int * l2g,
					 const double * grad_v0,
					 double * elementAverages,
					 int * tag,
					 double * Uin,
					 double * Uout);

extern
void applyDurlofskyDGlimiterP1Lagrange3d(int killExtrema,
					 int allowMinWithUndershoot,
					 int nElements_global,
					 int nElementBoundaries_element,
					 int nNodes_element,
					 int nSpace,
					 int nDOF_element,
					 const int * elementNeighborsArray,
					 const int * elementBoundariesArray,
					 const int * elementNodesArray,
					 const double * nodeArray,
					 const double * elementBarycentersArray,
					 const double * elementBoundaryBarycentersArray,
					 const double * elementNeighborShapeGradients,
					 const int * l2g,
					 const double * grad_v0,
					 double * elementAverages,
					 int * tag,
					 double * Uin,
					 double * Uout);
extern
void applyDurlofskyDGlimiterP1Lagrange2d_withVacuumTol(int killExtrema,
						       int allowMinWithUndershoot,
						       int enforcePositivity,
						       double vacuumTol,
						       int nElements_global,
						       int nElementBoundaries_element,
						       int nNodes_element,
						       int nSpace,
						       int nDOF_element,
						       const int * elementNeighborsArray,
						       const int * elementBoundariesArray,
						       const int * elementNodesArray,
						       const double * nodeArray,
						       const double * elementBarycentersArray,
						       const double * elementBoundaryBarycentersArray,
						       const double * elementNeighborShapeGradients,
						       const int * l2g,
						       const double * grad_v0,
						       double * elementAverages,
						       int * tag,
						       double * Uin,
						       double * Uout);

#endif
