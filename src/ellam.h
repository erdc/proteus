#ifndef ELLAM_H
#define ELLAM_H

#include <cassert>
#include <map>
#include <list>
#include <vector>
#include <iostream>
namespace ELLAM
{
inline
void evaluateTestFunctionsOnElement_linear_simplex_1d(const double x[3], const double x0[3], const double x1[3], double* w)
{
  double lambda_0 = (x[0]-x1[0])/(x0[0]-x1[0]);
  w[0] = lambda_0; w[1] = 1.0-lambda_0;
}
inline 
void evaluateTestFunctionsOnElement_linear_simplex_2d(const double x[3], 
						      const double x0[3], const double x1[3], const double x2[3],
						      const double n0[2], const double n1[2], const double n2[2],
						      double* w)
{
  //lambda_i = 1 - (x-a_i).n_i/(a_j-a_i).n_i
  double lambda_0 = 1.0 -
    ((x[0] - x0[0])*n0[0] + (x[1]-x0[1])*n0[1])/((x1[0]-x0[0])*n0[0] + (x1[1]-x0[1])*n0[1]);
  double lambda_1 = 1.0 -
    ((x[0] - x1[0])*n1[0] + (x[1]-x1[1])*n1[1])/((x0[0]-x1[0])*n1[0] + (x0[1]-x1[1])*n1[1]);

  w[0] = lambda_0; w[1] = lambda_1; w[2] = 1.0 - lambda_0 - lambda_1;
}
inline 
void evaluateTestFunctionsOnElement_linear_simplex_3d(const double x[3], 
						      const double x0[3], const double x1[3], const double x2[3], const double x3[3],
						      const double n0[3], const double n1[3], const double n2[3], const double n3[3],
						      double* w)
{
  //lambda_i = 1 - (x-a_i).n_i/(a_j-a_i).n_i
  double lambda_0 = 1.0 -
    ((x[0] - x0[0])*n0[0] + (x[1]-x0[1])*n0[1] + (x[2]-x0[2])*n0[2])/((x1[0]-x0[0])*n0[0] + (x1[1]-x0[1])*n0[1] + (x1[2]-x0[2])*n0[2]);
  double lambda_1 = 1.0 -
    ((x[0] - x1[0])*n1[0] + (x[1]-x1[1])*n1[1] + (x[2]-x1[2])*n1[2])/((x0[0]-x1[0])*n1[0] + (x0[1]-x1[1])*n1[1] + (x0[2]-x1[2])*n1[2]);
  double lambda_2 = 1.0 -
    ((x[0] - x2[0])*n2[0] + (x[1]-x2[1])*n2[1] + (x[2]-x2[2])*n2[2])/((x1[0]-x2[0])*n2[0] + (x1[1]-x2[1])*n2[1] + (x1[2]-x2[2])*n2[2]);
  
  w[0] = lambda_0; w[1] = lambda_1; w[2] = lambda_2; w[3] = 1.0 - lambda_0 - lambda_1 -lambda_2;
}

inline  
void evaluateTestFunctionsOnElement(int nSpace,
				    int nDOF_test_element,
				    int eN,
				    int nNodes_element,
				    int nElementBoundaries_element,
				    const double* nodeArray,
				    const int* elementNodesArray,
				    const double* elementBoundaryOuterNormalsArray,
				    const double x[3],
				    double* w)
{
  if (nSpace == 1)
    {
      assert(nDOF_test_element == 2); //linear only right now
      evaluateTestFunctionsOnElement_linear_simplex_1d(x,
						       &nodeArray[elementNodesArray[eN*nNodes_element+0]*3],
						       &nodeArray[elementNodesArray[eN*nNodes_element+1]*3],
						       w);
    }
  else if (nSpace == 2)
    {
      assert(nDOF_test_element == 3); //linear only right now
      evaluateTestFunctionsOnElement_linear_simplex_2d(x,
						       &nodeArray[elementNodesArray[eN*nNodes_element+0]*3],
						       &nodeArray[elementNodesArray[eN*nNodes_element+1]*3],
						       &nodeArray[elementNodesArray[eN*nNodes_element+2]*3],
						       &elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace+0*nSpace],
						       &elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace+1*nSpace],
						       &elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace+2*nSpace],
						       w);
    }
  else
    {
      assert(nDOF_test_element == 4); //linear only right now
      evaluateTestFunctionsOnElement_linear_simplex_3d(x,
						       &nodeArray[elementNodesArray[eN*nNodes_element+0]*3],
						       &nodeArray[elementNodesArray[eN*nNodes_element+1]*3],
						       &nodeArray[elementNodesArray[eN*nNodes_element+2]*3],
						       &nodeArray[elementNodesArray[eN*nNodes_element+3]*3],
						       &elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace+0*nSpace],
						       &elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace+1*nSpace],
						       &elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace+2*nSpace],
						       &elementBoundaryOuterNormalsArray[eN*nElementBoundaries_element*nSpace+3*nSpace],
						       w);

    }
}

inline
void exteriorOutflowFlux_c(int nSpace,
			   const double* n,
			   const double& u,
			   const double* velocity,
			   double& flux)
{
  flux = 0.0;
  double flow=0.0;
  for (int I=0; I < nSpace; I++)
    flow += n[I]*velocity[I];

  if (flow >= 0.0)
    {
      flux = u*flow;
    }
}

inline
void exteriorOutflowFluxDerivative_c(int nSpace,
				     const double* n,
				     const double* velocity,
				     double& dflux)
{
  double flow=0.0;
  for (int I=0; I < nSpace; I++)
    flow += n[I]*velocity[I];
  if (flow >= 0.0)
    dflux=flow;
  else
    dflux = 0.0;
}

inline
double ExteriorElementBoundaryFlux_c(const double& flux,
				     const double& w_dS)
{
  return flux*w_dS;
}

inline
double ExteriorNumericalAdvectiveFluxJacobian_c(const double& dflux_left,
						const double& v)
{
  return dflux_left*v;
}

//on output inflow >= -1 means track
//          
inline
void tagInflowPointForTracking_c(int nSpace,
				 const double& tn,
				 const double& tnp1,
				 const double& t,
				 const double* n,
				 const double* velocity_n,
				 const double* velocity_np1,
				 int& inflow)
{
  double flow=0.0;
  double dt = tnp1-tn;
  for (int I=0; I < nSpace; I++)
    flow += n[I]*((t-tn)/dt*velocity_np1[I] + (tnp1-t)/dt*velocity_n[I]);
  inflow = -3;
  if (flow < 0.0)
    {
      inflow = -1;
    }
}
double modifiedVanLeer(double r);
void collectTrackedPointsInElements(int nElements_global,
				   int nNodes_element,
				   int nElementBoundaries_element,
				   int nSpace,
				   int nPoints_tracked,
				   double boundaryTolerance,
				   double neighborTolerance,
				   const double* nodeArray,
				   const int* elementNodesArray,
				   const int* elementBoundariesArray,
				   const double* elementBoundaryLocalOuterNormalsArray,
				   const double* elementBoundaryBarycentersArray,
				   const int* element_track,
				   const int* flag_track,
				   const double* x_track,
				   std::map<int,std::list<double> >& elementsToTrackedPoints);


void createElementSSIPs_1d(const std::list<double>& elementTrackedPoints,
			   std::vector<double>& x_element, 
			   std::vector<double>& w_element);

};//namespace ELLAM



extern "C"
{
double updateOldMass_weak(int nSpace,            
			int nDOF_test_element, //dim for test function eval
			int nElements_global,  //mesh representation
			int nNodes_global,
			int nNodes_element,
			int nElementBoundaries_element,
			int nQuadraturePoints_element,     //element quadrature point data structures
			const double * nodeArray,          //mesh representation
			const int * elementNodesArray,
			const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
			const double * elementBoundaryOuterNormalsArray, //local element boundary outer normal constant on face
			const double* dV,             //integration weights 
			const double* x_track,        //location of forward tracked integration points, assumed for now 
			                                    //nElements_global x nQuadraturePoints_element
			const double* t_track,        //time forward tracked points stopped (t^{n+1} or earlier if exited domain)
			const int* element_track,     //element each forward tracked point ended up in
			const int* flag_track,        //id for each point, -1 -- interior, -2 exited domain, -3 didn't track for some reason
			const int* u_l2g,             //solution representation
			const double* q_m_last,
			double* q_elementResidual_u); 
  //int offset_u, int stride_u, 
  //		double* globalResidual);
void updateOldMass_weak_arbitraryQuadrature(int nSpace,            
					    int nDOF_test_element, //dim for test function eval
					    int nElements_global,  //mesh representation
					    int nNodes_global,
					    int nNodes_element,
					    int nElementBoundaries_element,
					    int nQuadraturePoints_track,     //element quadrature point data structures
					    const double * nodeArray,          //mesh representation
					    const int * elementNodesArray,
					    const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					    const double * elementBoundaryOuterNormalsArray, //local element boundary outer normal constant on face
					    const double* dV_track,             //integration weights at tracked points
					    const double* x_track,        //location of forward tracked integration points
					    const double* t_track,        //time forward tracked points stopped (t^{n+1} or earlier if exited domain)
					    const int* element_track,     //element each forward tracked point ended up in
					    const int* flag_track,        //id for each point, -1 -- interior, -2 exited domain, -3 didn't track for some reason
					    const int* u_l2g,             //solution representation
					    const double* q_m_track,     //mass from old time level evaluated at tracked points
					    double* q_elementResidual_u); 

double updateNewMass_weak(int nSpace,            
			  int nDOF_test_element, //dim for test function eval
			  int nElements_global,  //mesh representation
			  int nNodes_global,
			  int nNodes_element,
			  int nElementBoundaries_element,
			  int nQuadraturePoints_element,     //element quadrature point data structures
			  const double * nodeArray,          //mesh representation
			  const int * elementNodesArray,
			  const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
			  const double * elementBoundaryOuterNormalsArray, //local element boundary outer normal constant on face
			  const double* dV,             //integration weights 
			  const double* x,        //location of integration points, 
			                                    //nElements_global x nQuadraturePoints_element
			  const int* u_l2g,             //solution representation
			  const double* q_m,
			  double* q_elementResidual_u);

void evaluateSolutionAtTrackedPoints(int nSpace,
				     int nDOF_trial_element,
				     int nPoints_tracked,
				     int nElements_global,  //mesh representation
				     int nNodes_global,
				     int nNodes_element,
				     int nElementBoundaries_element,
				     const double * nodeArray,          //mesh representation
				     const int * elementNodesArray,
				     const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
				     const double * elementBoundaryOuterNormalsArray, //local element boundary outer normal constant on face
				     const double* x_track,        //location of  tracked integration points, assumed for now 
				                                   //nElements_global x nQuadraturePoints_element
				     const double* t_track,        //time tracked points stopped (t^{n+1} or earlier if exited domain)
       				     const int* element_track,     //element each tracked point ended up in
				     const int* flag_track,        //id for each point, -1 -- interior, -2 exited domain, -3 didn't track for some reason
				     const int* u_l2g,             //solution representation
				     const double* u_dof,           
				     double *u_x_track);            //u(x)


double updateExteriorOutflowBoundaryFlux(double dtnp1,          //full time step size
				       int nSpace,            
				       int nDOF_test_element, //dim for test function evalint nExteriorElementBoundaries_global,
				       int nQuadraturePoints_elementBoundary,
				       int nExteriorElementBoundaries_global,
				       const int* exteriorElementBoundariesArray,
				       const int* elementBoundaryElementsArray,
				       const int* elementBoundaryLocalElementBoundariesArray,
				       const double* ebqe_velocity_ext,
				       const double* ebqe_n_ext,
				       const double* ebqe_outflow_flux_last,
				       const double* u_test_dS_ext,
				       const double* ebqe_u,
				       const int* u_l2g,
				       double* ebqe_outflow_flux,
				       double* q_elementResidual_u); 

double updateExteriorOutflowBoundaryFluxInGlobalResidual(double dtnp1,          //full time step size
						       int nSpace,            
						       int nDOF_test_element, //dim for test function evalint nExteriorElementBoundaries_global,
						       int nQuadraturePoints_elementBoundary,
						       int nExteriorElementBoundaries_global,
						       const int* exteriorElementBoundariesArray,
						       const int* elementBoundaryElementsArray,
						       const int* elementBoundaryLocalElementBoundariesArray,
						       const double* ebqe_velocity_ext,
						       const double* ebqe_n_ext,
						       const double* ebqe_outflow_flux_last,
						       const double* u_test_dS_ext,
						       const double* ebqe_u,
						       const int* u_l2g,
						       double* ebqe_outflow_flux,
						       int offset_u, int stride_u, 
						       double* q_elementResidual_u, 
						       double* globalResidual);




void updateExteriorOutflowBoundaryFluxJacobian(double dtnp1,          //full time step size
					       int nSpace,            
					       int nDOF_test_element, //dim for test function evalint nExteriorElementBoundaries_global,
					       int nDOF_trial_element,
					       int nQuadraturePoints_elementBoundary,
					       int nExteriorElementBoundaries_global,
					       const int* exteriorElementBoundariesArray,
					       const int* elementBoundaryElementsArray,
					       const int* elementBoundaryLocalElementBoundariesArray,
					       const double* ebqe_velocity_ext,
					       const double* ebqe_n_ext,
					       const double* ebqe_outflow_flux_last,
					       const double* u_test_dS_ext,
					       const double* ebqe_u,
					       const double* u_trial_ext,
					       double* fluxJacobian);

void updateExteriorOutflowBoundaryFluxGlobalJacobian(double dtnp1,          //full time step size
						     int nSpace,            
						     int nDOF_test_element, //dim for test function evalint nExteriorElementBoundaries_global,
						     int nDOF_trial_element,
						     int nQuadraturePoints_elementBoundary,
						     int nExteriorElementBoundaries_global,
						     const int* exteriorElementBoundariesArray,
						     const int* elementBoundaryElementsArray,
						     const int* elementBoundaryLocalElementBoundariesArray,
						     const double* ebqe_velocity_ext,
						     const double* ebqe_n_ext,
						     const double* ebqe_outflow_flux_last,
						     const double* u_test_dS_ext,
						     const double* ebqe_u,
						     const double* u_trial_ext,
						     const int* csrRowIndeces_u_u, 
						     const int* csrColumnOffsets_u_u,
						     const int* csrColumnOffsets_eb_u_u,
						     double* globalJacobian);



void markInflowBoundaryPoints(int nSpace,
			      double tn, 
			      double tnp1,
			      double t,
			      int nExteriorElementBoundaries_global,
			      int nQuadraturePoints_elementBoundary,
			      const int* exteriorElementBoundariesArray,
			      const int* elementBoundaryElementsArray,
			      const int* elementBoundaryLocalElementBoundariesArray,
			      const double* ebqe_x,
			      const double* ebqe_n_ext,
			      const double* ebqe_velocity_ext_last,
			      const double* ebqe_velocity_ext,
			      const int* isDOFBoundary_u,
			      const int* isFluxBoundary_u,
			      int* element_track,//which element is point in
			      int* flag_track);  //>=-1 track


double accumulateInflowFlux(int nSpace,
			  int nDOF_test_element,
			  int nElements_global,//mesh representation
			  int nNodes_global,
			  int nNodes_element,
			  int nElementBoundaries_element,
			  int nExteriorElementBoundaries_global,
			  int nQuadraturePoints_elementBoundary,
			  const double * nodeArray,
			  const int * elementNodesArray,
			  const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
			  const int* exteriorElementBoundariesArray,
			  const int* elementBoundaryElementsArray,
			  const int* elementBoundaryLocalElementBoundariesArray,
			  const double * elementBoundaryOuterNormalsArray, //local element boundary outer normal constant on face
			  double tp, //time level for integration points
			  double timeWeight, //temporal integration weight (uniform for now)
			  const double* dS,  //spatial integration weights
			  const double* x_track_ext, //location of forward tracked integration points, assumed for now 
			                             //nExteriorElementBoundaries_global x nQuadraturePoints_elementBoundary
			  const double* t_track_ext,     //time forward tracked points stopped (t^{n+1} or earlier if exited domain)
			  const int* element_track_ext,  //element each forward tracked point ended up in
			  const int* flag_track_ext,     //id for each point, -1 -- interior, -2 exited domain, -3 didn't track for some reason
			  const int* u_l2g, 
			  const double* u_dof,
			  double* q_elementResidual_u, 
			  const int* sdInfo_u_rowptr, 
			  const int* sdInfo_u_colind,
			  const int* isFluxBoundary_u,
			  const double* ebqe_bc_flux_u_ext);

double accumulateInflowFluxInGlobalResidual(int nSpace,
					  int nDOF_test_element,
					  int nElements_global,//mesh representation
					  int nNodes_global,
					  int nNodes_element,
					  int nElementBoundaries_element,
					  int nExteriorElementBoundaries_global,
					  int nQuadraturePoints_elementBoundary,
					  const double * nodeArray,
					  const int * elementNodesArray,
					  const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					  const int* exteriorElementBoundariesArray,
					  const int* elementBoundaryElementsArray,
					  const int* elementBoundaryLocalElementBoundariesArray,
					  const double * elementBoundaryOuterNormalsArray, //local element boundary outer normal constant on face
					  double tp, //time level for integration points
					  double timeWeight, //temporal integration weight (uniform for now)
					  const double* dS,  //spatial integration weights
					  const double* x_track_ext, //location of forward tracked integration points, assumed for now 
					                             //nExteriorElementBoundaries_global x nQuadraturePoints_elementBoundary
					  const double* t_track_ext,     //time forward tracked points stopped (t^{n+1} or earlier if exited domain)
					  const int* element_track_ext,  //element each forward tracked point ended up in
					  const int* flag_track_ext,     //id for each point, -1 -- interior, -2 exited domain, -3 didn't track for some reason
					  const int* u_l2g, 
					  const double* u_dof,
					  double* q_elementResidual_u, 
					  const int* sdInfo_u_rowptr, 
					  const int* sdInfo_u_colind,
					  int offset_u, int stride_u, 
					  const int* isFluxBoundary_u,
					  const double* ebqe_bc_flux_u_ext,
					  double* globalResidual);


void tagNegligibleIntegrationPoints(int nPoints,
				    double zeroTol,
				    const double* x,
				    const double* u,
				    int* flag_track);


void calculateSlumpedMassApproximation1d(int nElements_global,
					 int nElementBoundaries_element,
					 int nQuadraturePoints_element,
					 int nDOF_trial_element,
					 int nDOF_test_element,
					 const int* l2g,
					 const int* elementNeighborsArray,
					 const double* u_dof,
					 const double* u_dof_limit,
					 const double* dm,
					 const double* w,
					 const double* v,
					 const double* dV,
					 const double* rhs,
					 double* elementResidual,
					 double* theta,
					 double* slumpedMassMatrixCorrection);
void calculateSlumpedMassApproximation1d_local(int nElements_global,
					       int nElementBoundaries_element,
					       int nQuadraturePoints_element,
					       int nDOF_trial_element,
					       int nDOF_test_element,
					       const int* l2g,
					       const int* elementNeighborsArray,
					       const double* u_dof,
					       const double* u_dof_limit,
					       const double* dm,
					       const double* w,
					       const double* v,
					       const double* dV,
					       const double* rhs,
					       double* elementResidual,
					       double* theta,
					       double* slumpedMassMatrixCorrection);

void calculateSlumpedMassApproximation2dOrig(int nElements_global,
					 int nElementBoundaries_element,
					 int nQuadraturePoints_element,
					 int nDOF_trial_element,
					 int nDOF_test_element,
					 const int* l2g,
					 const int* elementNeighborsArray,
					 const double* u_dof,
					 const double* u_dof_limit,
					 const double* dm,
					 const double* w,
					 const double* v,
					 const double* dV,
					 const double* rhs,
					 double* elementResidual,
					 double* theta,
					 double* slumpedMassMatrixCorrection);
void calculateSlumpedMassApproximation2d(int nElements_global,
					 int nElementBoundaries_element,
					 int nQuadraturePoints_element,
					 int nDOF_trial_element,
					 int nDOF_test_element,
					 double adjustFactor,
					 const int* l2g,
					 const int* elementNeighborsArray,
					 const double* u_dof,
					 const double* u_dof_limit,
					 const double* dm,
					 const double* w,
					 const double* v,
					 const double* dV,
					 const double* rhs,
					 double* elementResidual,
					 double* theta,
					 double* slumpedMassMatrixCorrection);

void calculateSlumpedMassApproximation2d_upwind(int nElements_global,
						int nNodes_global,
						int nElementBoundaries_element,
						int nNodes_element,
						int nQuadraturePoints_element,
						int nDOF_trial_element,
						int nDOF_test_element,
						const double* nodeArray,
						const int* elementNodesArray,
						const int* elementNeighborsArray,
						const int* nodeStarArray,
						const int* nodeStarOffsets,
						const double* elementBoundaryLocalOuterNormalsArray,
						const int* l2g,
						const double* u_dof,
						const double* u_dof_limit,
						const double* dm,
						const double* df,
						const double* w,
						const double* v,
						const double* dV,
						const double* rhs,
						double* elementResidual,
						double* theta,
						double* slumpedMassMatrixCorrection);


void updateElementJacobianWithSlumpedMassApproximation(int nElements_global,
						       int nDOF_trial_element,
						       int nDOF_test_element,
						       const double* theta,
						       double* elementJacobianMatrix);
void calculateBerzinsSlumpedMassApproximation1d(int nElements_global,
					       int nElementBoundaries_element,
					       int nQuadraturePoints_element,
					       int nDOF_trial_element,
					       int nDOF_test_element,
					       const int* l2g,
					       const int* elementNeighborsArray,
					       const double* u_dof,
					       const double* u_dof_limit,
					       const double* dm,
					       const double* w,
					       const double* v,
					       const double* dV,
					       const double* rhs,
					       double* elementResidual,
					       double* slumpedMassMatrixCorrection);
void calculateBerzinsSlumpedMassApproximation1dv2(int nElements_global,
					       int nElementBoundaries_element,
					       int nQuadraturePoints_element,
					       int nDOF_trial_element,
					       int nDOF_test_element,
					       const int* l2g,
					       const int* elementNeighborsArray,
					       const double* u_dof,
					       const double* u_dof_limit,
					       const double* dm,
					       const double* w,
					       const double* v,
					       const double* dV,
					       const double* rhs,
					       double* elementResidual,
					       double* slumpedMassMatrixCorrection);
void calculateBerzinsSlumpedMassApproximation2d(int nElements_global,
						int nElementBoundaries_element,
						int nQuadraturePoints_element,
						int nDOF_trial_element,
						int nDOF_test_element,
						const int* l2g,
						const int* elementNeighborsArray,
						const double* u_dof,
						const double* u_dof_limit,
						const double* dm,
						const double* w,
						const double* v,
						const double* dV,
						const double* rhs,
						double* elementResidual,
						double* slumpedMassMatrixCorrection);
void updateElementJacobianWithSlumpedMassCorrection(int nElements_global,
						    int nDOF_trial_element,
						    int nDOF_test_element,
						    const double* elementMassMatrixCorrection,
						    double* elementJacobianMatrix);


void manuallyUpdateGlobalMassMatrix(int nElements_global,
				    int nQuadraturePoints_element,
				    int nDOF_trial_element,
				    int nDOF_test_element,
				    const int* rowptr,
				    const int* colind,
				    const int* l2g,
				    const double* u_dof,
				    const double* dm,
				    const double* w,
				    const double* v,
				    const double* dV,
				    double* globalMassMatrix);

void calculateElementSlumpedMassApproximationFromGlobalEdgeLimiter(int nElements_global,
								   int nQuadraturePoints_element,
								   int nDOF_trial_element,
								   int nDOF_test_element,
								   const int* rowptr,
								   const int* colind,
								   const int* l2g,
								   const double* u_dof,
								   const double* dm,
								   const double* w,
								   const double* v,
								   const double* dV,
								   const double* globalEdgeSlumpingParameter,
								   double* elementResidual,
								   double* slumpedMassMatrixCorrection);

void computeSlumpingParametersFCT_KuzminMoeller10(const int nDOF_global,
						  const int* rowptr, //sparsity information
						  const int* colind,
						  const double* u_dof,
						  const double* u_dof_limit,
						  const double* Mc, //assembled global mass matrix
						  double * Rip,  //waste some temporary space until get 
						  double * Rim,  //algorithm debugged
						  double* edgeSlumpingParameter); //spare rep for slumping factor for edge ij (i.e. alpha_ij)

  //SSIP stuf
void generateQuadratureArraysForSSIPs(int nElements_global,
			    int nNodes_element,
			    int nElementBoundaries_element,
			    int nSpace,
			    int nPoints_tracked,
			    int nPointsPerElement_default,
			    double boundaryTolerance,
			    double neighborTolerance,
			    const double* nodeArray,
			    const int* elementNodesArray,
			    const int* elementBoundariesArray,
			    const double* elementBoundaryLocalOuterNormalsArray,
			    const double* elementBoundaryBarycentersArray,
			    const int* element_track,
			    const int* flag_track,
			    const double* x_track,
			    const double* q_x_default,
			    const double* q_dV_default,
			    std::vector<int>& elements_gq,
			    std::vector<double>& x_qg,
			    std::vector<double>& dV_gq);

void generateArraysForTrackedSSIPs(int nElements_global,
				   int nNodes_element,
				   int nElementBoundaries_element,
				   int nSpace,
				   int nPoints_tracked,
				   double boundaryTolerance,
				   double neighborTolerance,
				   const double* nodeArray,
				   const int* elementNodesArray,
				   const int* elementBoundariesArray,
				   const double* elementBoundaryLocalOuterNormalsArray,
				   const double* elementBoundaryBarycentersArray,
				   const int* element_track,
				   const int* flag_track,
				   const double* x_track,
				   std::vector<int>& element_offsets_ssip,
				   std::vector<double>& x_ssip);


double integratePiecewiseLinearMassSource(int nknots,
					 const double * t_vals,
					 const double * m_vals,
					 double t_in,
					 double t_out,
					 double tau_out,
					 double tau_in,
					 double decay); //must be nonpositive

void accumulateSourceContribution(int nParticles_global, //number of particles in this source
				  int nElements_global,  //number of elements in domain
				  int nParticleFlags,    //total number of particle types or flags
				  int nMassSourceKnots,  //number of knots in source spline
				  double tau,            //time evaluating solution
				  const int * traj_offsets, //traj_offsets[i] = start of trajectory info for particle i
				                            //n_i = traj_offsets[i+1]-traj_offsets[i] 
				  const double * x_traj, //particle trajectories: x (size 3\sum_i n_i])
				  const double * t_traj, //particle trajectories: t
				  const int * elem_traj, //particle trajectories: element id's
				  const double * massSource_t, //discrete t values (knot's) for mass source
				  const double * massSource_m, //values for mass source at knot's
				  const double * decay_coef_element,  //linear decay: nParticleFlags * nElements_global
				  const double * retardation_factor_element,  //Retardation factor: nParticleFlags * nElements_global
				  const int * particleFlags, //The particle 'type' associated with particles
				  double *c_element); //element concentrations

void accumulateSourceContributionMaterialTypes(int nParticles_global, //number of particles in this source
					       int nElements_global,  //number of elements in domain
					       int nParticleFlags,    //total number of particle types or flags
					       int nMassSourceKnots,  //number of knots in source spline
					       double tau,            //time evaluating solution
					       const int * traj_offsets, //traj_offsets[i] = start of trajectory info for particle i
  				                                         //n_i = traj_offsets[i+1]-traj_offsets[i]  
					       const double * x_traj, //particle trajectories: x (size 3\sum_i n_i])
					       const double * t_traj, //particle trajectories: t
					       const int * elem_traj, //particle trajectories: element id's
					       const double * massSource_t, //discrete t values (knot's) for mass source
					       const double * massSource_m, //values for mass source at knot's
					       const int* material_types_element, //identifier for material of each element 
					       const double * decay_coef_types,  //linear decay:  nMaterialTypes * nParticleFlags 
					       const double * retardation_factor_types,  //Retardation factor: nMaterialTypes * nParticleFlags  
					       const int * particle_flags, //The particle 'type' associated with particles
					       double *c_element); //element concentrations

double volume123(int nSpace, //space dimension
		 int nNodes_element, //number of nodes in this element
		 const int* elementNodes_element, //element--> global node table for this element 
		 const double* nodeArray);

void nodalProjection123(int calculateNodeStar_volume,
			int nElements_global,
			int nNodes_global,
			const int* elementNodes_offset,
			const int* elementNodesArray,
			const double* elementVolumesArray,
			const double* c_ele,
			double * nodeStar_volume, 
			double * c_nodal);
}//extern C
#endif 
