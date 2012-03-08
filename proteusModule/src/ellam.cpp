#include "ellam.h"
#include "tracking.h"
#include <cmath>
#include <iostream>
#include <cassert>
#include <list>
#include <map>
#include <algorithm>
#include <string.h>

#include PROTEUS_LAPACK_H
/***********************************************************************

  TODO:

   switch dimension of dV and x_track to just be the number of global
   integration points so that there can be different numbers per element

   decide if ignoring outflow mass is a better idea
**********************************************************************/

extern "C" 
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
			  double* q_elementResidual_u) 
{
  using namespace ELLAM;
  /***********************************************************************
    for current integration point x_k
       1. get mass at old time level, m^{n}_k
       2. get integration weight W_k
       3. get forward tracked integration point x^{n+1}_k
       4. evaluate test functions with non-zero support at x^{n+1}_k
           a. get element eN^{n+1}_k from element_track[k]
           b. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       5. accumulate m^{n+1}_k*W_k*w_i(x^{n+1}_{k}) to element residual 
  **********************************************************************/
  /*return value*/
  double totalOldMass = 0.0;
  //for now just assume a max dim
  register double w[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  assert(nDOF_test_element <= 10);
  for(int eN=0;eN<nElements_global;eN++)
    {
      //loop over quadrature points, may end up on elements with different support
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      w[i] = 0.0;
	    }
	  //compute indeces and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k;
	  //todo decide if ignoring outflow mass is a better idea
	  if (flag_track[eN_k] >= -2) //-2 exited domain
	    {
	      //m_k^{n+1},W_k
	      register double m_old = q_m_last[eN_k],
		weight = dV[eN_k];
	      register int eN_track = element_track[eN_k];

	      //mwf debug
// 	      std::cout<<"ellam tracked point ["<<eN<<","<<k<<"] --> "<<x_track[eN_k*3]
// 	      	       <<" eN_track= "<<eN_track
// 		       <<" x0= ["<<nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3]
// 		       <<","<<nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3+1]<<"]"
// 		       <<" x1= ["<<nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3]
// 		       <<","<<nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3+1]<<"]"
// 		       <<" x2= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+2]*3]
// 		       <<","<<nodeArray[elementNodesArray[eN_track*nNodes_element+2]*3+1]<<"]"<<std::endl
// 		       <<" n0= ["<<elementBoundaryOuterNormalsArray[eN_track*nElementBoundaries_element*nSpace+0*nSpace]
// 		       <<","<<elementBoundaryOuterNormalsArray[eN_track*nElementBoundaries_element*nSpace+0*nSpace+1]<<"]"
// 		       <<" n1= ["<<elementBoundaryOuterNormalsArray[eN_track*nElementBoundaries_element*nSpace+1*nSpace]
// 		       <<","<<elementBoundaryOuterNormalsArray[eN_track*nElementBoundaries_element*nSpace+1*nSpace+1]<<"]"
// 		       <<" n2= ["<<elementBoundaryOuterNormalsArray[eN_track*nElementBoundaries_element*nSpace+2*nSpace]
// 		       <<","<<elementBoundaryOuterNormalsArray[eN_track*nElementBoundaries_element*nSpace+2*nSpace+1]<<"]"
// 		       <<std::endl;
	      evaluateTestFunctionsOnElement(nSpace,
					     nDOF_test_element,
					     eN_track,
					     nNodes_element,
					     nElementBoundaries_element,
					     nodeArray,
					     elementNodesArray,
					     elementBoundaryOuterNormalsArray,
					     &x_track[eN_k*3],
					     w);
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int eN_track_i=eN_track*nDOF_test_element+i;
		  q_elementResidual_u[eN_track_i] -= m_old*w[i]*weight;
		  totalOldMass += m_old*w[i]*weight;
		  //globalResidual[offset_u+stride_u*u_l2g[eN_track_i]] -= m_old*w[i]*weight; 
		  //mwf debug
		  //std::cout<<"LADRellam oldMass globalResid["<<offset_u+stride_u*u_l2g[eN_track_i]<<"]= "<< globalResidual[offset_u+stride_u*u_l2g[eN_track_i]] 
		  //   <<" oldmass = "<<m_old*w[i]*weight <<std::endl;
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element
    }//nElements loop
  return totalOldMass;
}//updateOldMass
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
					    double* q_elementResidual_u) 
{
  using namespace ELLAM;
  /***********************************************************************
    for current integration point x_k
       1. get mass at old time level, m^{n}_k
       2. get integration weight W_k
       3. get forward tracked integration point x^{n+1}_k
       4. evaluate test functions with non-zero support at x^{n+1}_k
           a. get element eN^{n+1}_k from element_track[k]
           b. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       5. accumulate m^{n+1}_k*W_k*w_i(x^{n+1}_{k}) to element residual 
  **********************************************************************/
  //for now just assume a max dim
  register double w[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  assert(nDOF_test_element <= 10);
  for(int k=0; k < nQuadraturePoints_track; k++)
    {
      for (int i=0;i<nDOF_test_element;i++)
	{
	  w[i] = 0.0;
	}
      //todo decide if ignoring outflow mass is a better idea
      if (flag_track[k] >= -1) //-2 is exited domain
	{
	  register int eN_track = element_track[k];
	  //m_k^{n+1},W_k
	  register double m_old = q_m_track[k],
	    weight = dV_track[k];

	  evaluateTestFunctionsOnElement(nSpace,
					 nDOF_test_element,
					 eN_track,
					 nNodes_element,
					 nElementBoundaries_element,
					 nodeArray,
					 elementNodesArray,
					 elementBoundaryOuterNormalsArray,
					 &x_track[k*3],
					 w);
	      
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_track_i=eN_track*nDOF_test_element+i;
	      q_elementResidual_u[eN_track_i] -= m_old*w[i]*weight;
	    }//i
	      
	}//needed to track point
	    
    }//integration points
}//updateOldMass_weak_arbitraryQuadrature

//include a function for evaluating new mass to see how much impact having different test function evaluates for
//new and old makes
extern "C" 
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
			  double* q_elementResidual_u) 
{
  using namespace ELLAM;
  /***********************************************************************
    for current integration point x_k on element eN
       1. get mass at new time level, m^{n+1}_k
       2. get integration weight W_k
       3. get integration point x^{n+1}_k
       4. evaluate test functions with non-zero support at x^{n+1}_k
           a. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       5. accumulate m^{n+1}_k*W_k*w_i(x^{n+1}_{k}) to element residual 
  **********************************************************************/
  /*return value*/
  double totalNewMass = 0.0;
  //for now just assume a max dim
  register double w[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  assert(nDOF_test_element <= 10);
  for(int eN=0;eN<nElements_global;eN++)
    {
      //loop over quadrature points, may end up on elements with different support
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      w[i] = 0.0;
	    }
	  //compute indeces and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k;
	  //m_k^{n+1},W_k
	  register double m_new = q_m[eN_k],
	    weight = dV[eN_k];
	  
	  evaluateTestFunctionsOnElement(nSpace,
					 nDOF_test_element,
					 eN,
					 nNodes_element,
					 nElementBoundaries_element,
					 nodeArray,
					 elementNodesArray,
					 elementBoundaryOuterNormalsArray,
					 &x[eN_k*3],
					 w);
	  
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_i=eN*nDOF_test_element+i;
	      q_elementResidual_u[eN_i] += m_new*w[i]*weight;
	      totalNewMass += m_new*w[i]*weight;
	    }//i
	}//integration point per element
    }//nElements loop
  return totalNewMass;
}//updateOldMass

extern "C" 
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
				     double *u_x_track)            //u(x)
{
  using namespace ELLAM;
  /***********************************************************************
    Evaluate solution at tracked points using physical space representation
      of solution for simplicity

    for current integration point x_k
       3. evaluate trial functions with non-zero support at x_k
           a. get element eN^{n}_k from element_track[k]
           b. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       4. get solution, u(x_k) = \sum u^j w_j(x_k)
  **********************************************************************/
  //for now just assume a max dim
  register double w[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  assert(nDOF_trial_element <= 10);
  for (int k=0; k < nPoints_tracked; k++)
    {
      register double u=0;
      for (int i=0;i<nDOF_trial_element;i++)
	{
	  w[i] = 0.0;
	}
      
      if (flag_track[k] == -1) //do not evaluate points that go out of domain, let inflow boundary approx take care of this?
	{
	  register int eN_track = element_track[k];
	  assert(eN_track >= 0 && eN_track < nElements_global);

	  //same as trial
	  evaluateTestFunctionsOnElement(nSpace,
					 nDOF_trial_element,
					 eN_track,
					 nNodes_element,
					 nElementBoundaries_element,
					 nodeArray,
					 elementNodesArray,
					 elementBoundaryOuterNormalsArray,
					 &x_track[k*3],
					 w);
	      
	  for(int i=0;i<nDOF_trial_element;i++) 
	    { 
	      u += w[i]*u_dof[u_l2g[eN_track*nDOF_trial_element+i]];
	    }//i
	  
	}//needed to track point
      u_x_track[k] = u;
    }//points loop
}//evaluateSolutionAtTrackedPoints

extern "C"
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
					 double* q_elementResidual_u)
{
  using namespace ELLAM;
  double totalOutflow = 0.0;
  /***************************************************
     apply outflow flux with simple trapezoidal rule
     in time
   **************************************************/

  //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
  //
  //ebNE is the Exterior element boundary INdex
  //ebN is the element boundary INdex
  //eN is the element index
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      register int ebN = exteriorElementBoundariesArray[ebNE], 
	eN  = elementBoundaryElementsArray[ebN*2+0];
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  register double u_ext=0.0,
	    flux_ext=0.0;
	  
	  u_ext =  ebqe_u[ebNE_kb];
	  
	  // 
	  //calculate outflow numerical fluxes, 
	  //for now applies zero diffusive flux if v.n > 0 
	  // 
	  exteriorOutflowFlux_c(nSpace,
				&ebqe_n_ext[ebNE_kb_nSpace],
				u_ext,
				&ebqe_velocity_ext[ebNE_kb_nSpace],
				flux_ext);
	  //save for next time step
	  ebqe_outflow_flux[ebNE_kb] = flux_ext;
	  //add time scaling
	  flux_ext = dtnp1*0.5*(flux_ext + ebqe_outflow_flux_last[ebNE_kb]);

	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i,
		eN_i = eN*nDOF_test_element+i;

	      q_elementResidual_u[eN_i] += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
	      totalOutflow += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
	      //mwf debug
	      //std::cout<<"LADR2Dellam boundaryResidual globalResid["<<offset_u+stride_u*u_l2g[eN_i]<<"]= "<< globalResidual[offset_u+stride_u*u_l2g[eN_i]] 
	      //   <<" elementResidual = "<<elementResidual_u[i] <<std::endl;
	    }//i
	}//kb
    }//ebNE
  return totalOutflow;
}
extern "C"
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
							 double* globalResidual)
{
  using namespace ELLAM;
  double totalOutflowFlux  = 0.0; /*return value*/
  /***************************************************
     apply outflow flux with simple trapezoidal rule
     in time
   **************************************************/

  //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
  //
  //ebNE is the Exterior element boundary INdex
  //ebN is the element boundary INdex
  //eN is the element index
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      register int ebN = exteriorElementBoundariesArray[ebNE], 
	eN  = elementBoundaryElementsArray[ebN*2+0];
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  register double u_ext=0.0,
	    flux_ext=0.0;
	  
	  u_ext =  ebqe_u[ebNE_kb];
	  
	  // 
	  //calculate outflow numerical fluxes, 
	  //for now applies zero diffusive flux if v.n >= 0 
	  // 
	  exteriorOutflowFlux_c(nSpace,
				&ebqe_n_ext[ebNE_kb_nSpace],
				u_ext,
				&ebqe_velocity_ext[ebNE_kb_nSpace],
				flux_ext);
	  //save for next time step
	  ebqe_outflow_flux[ebNE_kb] = flux_ext;
	  //add time scaling
	  flux_ext = dtnp1*0.5*(flux_ext + ebqe_outflow_flux_last[ebNE_kb]);

	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i,
		eN_i = eN*nDOF_test_element+i;

	      q_elementResidual_u[eN_i] += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
	      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
	      totalOutflowFlux += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
	      //mwf debug
	      //std::cout<<"LADR2Dellam boundaryResidual globalResid["<<offset_u+stride_u*u_l2g[eN_i]<<"]= "<< globalResidual[offset_u+stride_u*u_l2g[eN_i]] 
	      //   <<" elementResidual = "<<elementResidual_u[i] <<std::endl;
	    }//i
	}//kb
    }//ebNE
  return totalOutflowFlux;
}

extern "C"
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
						     double* globalJacobian)
{
  using namespace ELLAM;
  const int nDOF_test_X_trial_element = nDOF_test_element*nDOF_trial_element;
  /***************************************************
     apply outflow flux with simple trapezoidal rule
     in time
   **************************************************/
  //
  //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
  //
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      register int ebN = exteriorElementBoundariesArray[ebNE]; 
      register int eN  = elementBoundaryElementsArray[ebN*2+0];
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  register double u_ext=0.0,
	    dflux_u_u_ext=0.0;
	  
	  u_ext =  ebqe_u[ebNE_kb];

	  // 
	  //calculate the numerical fluxes 
	  // 
	  exteriorOutflowFluxDerivative_c(nSpace,
					  &ebqe_n_ext[ebNE_kb_nSpace],
					  &ebqe_velocity_ext[ebNE_kb_nSpace],
					  dflux_u_u_ext);

	  //
	  //make sure have correct dt scaling for boundaries
	  //trapezoidal rule so use 0.5dt
	  dflux_u_u_ext *= 0.5*dtnp1;
	  //
	  //update the global Jacobian from the flux Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i = eN*nDOF_test_element+i,
		ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
		  register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
		  register double fluxJacobian_u_u_j = ExteriorNumericalAdvectiveFluxJacobian_c(dflux_u_u_ext,u_trial_ext[ebNE_kb_j]);
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u_j*u_test_dS_ext[ebNE_kb_i];
		  //mwf debug
		  //std::cout<<"LADR2Dellam eN= "<<eN<<" ebN= "<<ebN<<" fluxJacobian_u_u["<<j<<"]= "<<fluxJacobian_u_u[j]<<std::endl;
		}//j
	    }//i
	}//kb
    }//ebNE
}//outflowFluxJacobian
extern "C"
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
					       double* fluxJacobian)
{
  using namespace ELLAM;
  /***************************************************
     apply outflow flux with simple trapezoidal rule
     in time
   **************************************************/
  //
  //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
  //
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  register double u_ext=0.0,
	    dflux_u_u_ext=0.0;
	  
	  u_ext =  ebqe_u[ebNE_kb];

	  // 
	  //calculate the numerical fluxes 
	  // 
	  exteriorOutflowFluxDerivative_c(nSpace,
					  &ebqe_n_ext[ebNE_kb_nSpace],
					  &ebqe_velocity_ext[ebNE_kb_nSpace],
					  dflux_u_u_ext);

	  //
	  //make sure have correct dt scaling for boundaries
	  //trapezoidal rule so use 0.5dt
	  dflux_u_u_ext *= 0.5*dtnp1;
	  //
	  //update the local flux Jacobian
	  //
	  for (int j=0; j < nDOF_trial_element; j++)
	    {
	      register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
	      register double fluxJacobian_u_u_j = ExteriorNumericalAdvectiveFluxJacobian_c(dflux_u_u_ext,u_trial_ext[ebNE_kb_j]);
	      fluxJacobian[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+
			   kb*nDOF_trial_element+
			   j]
		+= fluxJacobian_u_u_j;
	    }
	}//kb
    }//ebNE
}//outflowFluxJacobian

/**********************************************************************
  just tag points that need to be integrated forward for inflow boundary
  need to decide how to handle transient velocity field. Probably
  pass in new time level and old time level one next 
  Right now ignores DOF or Flux boundary flag, selects for tracking if
  velocity is inflow
 **********************************************************************/
extern "C" void markInflowBoundaryPoints(int nSpace,
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
					 int* flag_track)  //>=-1 track, -1 don't

{
  using namespace ELLAM;
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      const int ebN = exteriorElementBoundariesArray[ebNE];
      const int eN  = elementBoundaryElementsArray[ebN*2+0];
      for (int kb=0; kb < nQuadraturePoints_elementBoundary; kb++)
	{
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  int needToTrackPoint=0;
	  tagInflowPointForTracking_c(nSpace,
				      tn,tnp1,t,
				      &ebqe_n_ext[ebNE_kb_nSpace],
				      &ebqe_velocity_ext_last[ebNE_kb_nSpace],
				      &ebqe_velocity_ext[ebNE_kb_nSpace],
				      needToTrackPoint);
	  flag_track[ebNE_kb] = needToTrackPoint;
	  element_track[ebNE_kb]= eN;
	}
    }

}

/***********************************************************************
   set flag to -2 for points that have |u| < tol
 ***********************************************************************/
extern "C" void tagNegligibleIntegrationPoints(int nPoints,
					       double zeroTol,
					       const double* x,
					       const double* u,
					       int *flag_track)
{
  for (int k=0; k < nPoints; k++)
    {
      if (fabs(u[k]) < zeroTol)
	flag_track[k] = -2;
    }
}
					       
/**********************************************************************
  go through inflow boundary points that have been tracked forward
  and accumulate the inflow flux into the residual for the test functions
  with non-zero support in those areas
 **********************************************************************/
extern "C" 
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
					  double* globalResidual)
{

  using namespace ELLAM;
  /***********************************************************************
    for current boundary integration point (x_k,t_p)
       1. get boundary flux, \sigma^{p}_k
       2. get spatial integration weight W_k, temporal weight \Delta t^p
       3. get forward tracked integration point x^{n+1}_k
       4. evaluate test functions with non-zero support at x^{n+1}_k
           a. get element eN^{n+1}_k from element_track_ext[k]
           b. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       5. accumulate \sigma^{p}_k*\Delta t^p* W_k*w_i(x^{n+1}_{k}) to global residual I(i) 
  **********************************************************************/
  double totalInflow = 0.0;
  double w[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  assert(nDOF_test_element <= 10);
  //todo switch to a loop over boundary integration points?
  for(int ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      register int ebN = exteriorElementBoundariesArray[ebNE];
      
      //loop over quadrature points, may end up on elements with different support
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
        {
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      w[i] = 0.0;
	    }
	  //compute indeces and declare local storage
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb;
	  if (flag_track_ext[ebNE_kb] >= -1 && isFluxBoundary_u[ebNE_kb])
	    {
	      //sigma_k,W_k
	      register double totalFlux = ebqe_bc_flux_u_ext[ebNE_kb],
		spatialWeight = dS[ebNE_kb];
	      register int eN_track = element_track_ext[ebNE_kb];
	      //mwf debug
	      //std::cout<<"ellam tracked boundary point ["<<ebNE<<","<<kb<<"] --> "<<x_track[ebNE_kb*3]
	      //	       <<" eN_track= "<<eN_track<<" x0= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3]
	      //       <<" x1= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3]<<std::endl;
	      evaluateTestFunctionsOnElement(nSpace,
					     nDOF_test_element,
					     eN_track,
					     nNodes_element,
					     nElementBoundaries_element,
					     nodeArray,
					     elementNodesArray,
					     elementBoundaryOuterNormalsArray,
					     &x_track_ext[ebNE_kb*3],
					     w);
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int eN_track_i=eN_track*nDOF_test_element+i;
		  q_elementResidual_u[eN_track_i] += totalFlux*w[i]*spatialWeight*timeWeight;
		  globalResidual[offset_u+stride_u*u_l2g[eN_track_i]] += totalFlux*w[i]*spatialWeight*timeWeight; 
		  totalInflow += totalFlux*w[i]*spatialWeight*timeWeight;
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element boundary
    }//element boundaries
  return totalInflow;
}


/**********************************************************************
  go through inflow boundary points that have been tracked forward
  and accumulate the inflow flux into the element residual for the test functions
  with non-zero support in those areas
 **********************************************************************/
extern "C" 
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
			  const double* ebqe_bc_flux_u_ext)
{

  using namespace ELLAM;
  /***********************************************************************
    for current boundary integration point (x_k,t_p)
       1. get boundary flux, \sigma^{p}_k
       2. get spatial integration weight W_k, temporal weight \Delta t^p
       3. get forward tracked integration point x^{n+1}_k
       4. evaluate test functions with non-zero support at x^{n+1}_k
           a. get element eN^{n+1}_k from element_track_ext[k]
           b. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       5. accumulate \sigma^{p}_k*\Delta t^p* W_k*w_i(x^{n+1}_{k}) to global residual I(i) 
  **********************************************************************/
  double totalInflowFlux = 0.0;

  double w[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  assert(nDOF_test_element <= 10);
  //todo switch to a loop over boundary integration points?
  for(int ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      register int ebN = exteriorElementBoundariesArray[ebNE];
      
      //loop over quadrature points, may end up on elements with different support
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
        {
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      w[i] = 0.0;
	    }
	  //compute indeces and declare local storage
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb;
	  if (flag_track_ext[ebNE_kb] >= -1 && isFluxBoundary_u[ebNE_kb])
	    {
	      //sigma_k,W_k
	      register double totalFlux = ebqe_bc_flux_u_ext[ebNE_kb],
		spatialWeight = dS[ebNE_kb];
	      register int eN_track = element_track_ext[ebNE_kb];
	      //mwf debug
	      //std::cout<<"ellam tracked boundary point ["<<ebNE<<","<<kb<<"] --> "<<x_track[ebNE_kb*3]
	      //	       <<" eN_track= "<<eN_track<<" x0= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3]
	      //       <<" x1= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3]<<std::endl;
	      evaluateTestFunctionsOnElement(nSpace,
					     nDOF_test_element,
					     eN_track,
					     nNodes_element,
					     nElementBoundaries_element,
					     nodeArray,
					     elementNodesArray,
					     elementBoundaryOuterNormalsArray,
					     &x_track_ext[ebNE_kb*3],
					     w);
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int eN_track_i=eN_track*nDOF_test_element+i;
		  q_elementResidual_u[eN_track_i] += totalFlux*w[i]*spatialWeight*timeWeight;
		  totalInflowFlux += totalFlux*w[i]*spatialWeight*timeWeight;
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element boundary
    }//element boundaries
  return totalInflowFlux;
}


extern "C"
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
					 double* slumpedMassMatrixCorrection)
{
  int eN,i,j,k,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int I0,Im1,Ip1,Ip2;
  double drm1,dr0,drp1,a,b,tmp1,tmp2,tmp3,tmp4;
  int eN_neighbor0,eN_neighbor1;
  assert(nDOF_test_element == 2);
  for (eN=0; eN<nElements_global; eN++)
    {
      I0 = l2g[eN*nDOF_test_element + 0];
      Ip1= l2g[eN*nDOF_test_element + 1];
      dr0 = fabs(rhs[Ip1]-rhs[I0]);
      
      /*local stencil relative to node I*/
      /*I0 neighbor across from node I1 and vice versa*/
      eN_neighbor0 =elementNeighborsArray[eN*nElementBoundaries_element+1];
      eN_neighbor1 =elementNeighborsArray[eN*nElementBoundaries_element+0];

      
      drm1 = 0.0; Im1 = -1;
      if (eN_neighbor0 >= 0)
	{
	  for (j=0; j < nDOF_test_element; j++)
	    {
	      if (l2g[eN_neighbor0*nDOF_test_element + j] != I0)
		{
		  Im1 = l2g[eN_neighbor0*nDOF_test_element + j];
		  drm1 = fabs(rhs[I0]-rhs[Im1]);
		  break;
		}
	    }
	}
	  
      drp1 = 0.0; Ip2 = -1;
      if (eN_neighbor1 >= 0)
	{
	  for (j=0; j < nDOF_test_element; j++)
	    {
	      if (l2g[eN_neighbor1*nDOF_test_element + j] != Ip1)
		{
		  Ip2 = l2g[eN_neighbor1*nDOF_test_element + j];
		  drp1 = fabs(rhs[Ip2]-rhs[Ip1]);
		  break;
		}
	    }
	  
	}
      
      /*calculate local mass matrix entries isn't exactly right
	for now, since assumes matrix same on both elements*/
      a = 0.0; b = 0.0; i = 0; j = 1;
      for (k=0; k < nQuadraturePoints_element; k++)
	{
	  a += 
	    dm[eN*nQuadraturePoints_element+k]
	    *
	    v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
	      i]
	    *
	    w[eN*nQuadraturePoints_element*nDOF_trial_element + 
	      k*nDOF_test_element + 
	      i]
	    *
	    dV[eN*nQuadraturePoints_element+k];
	  b += 
	    dm[eN*nQuadraturePoints_element+k]
	    *
	    v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
	      j]
	    *
	    w[eN*nQuadraturePoints_element*nDOF_trial_element + 
	      k*nDOF_test_element + 
	      i]
	    *
	    dV[eN*nQuadraturePoints_element+k];
	}
      
      tmp1 = dr0/(drm1+dr0+drp1+1.0e-12);
      tmp2 = fmin(drm1,dr0)/(drm1+dr0+1.0e-12);
      tmp3 = fmin(drp1,dr0)/(drp1+dr0+1.0e-12);
      tmp4 = fmin(fmin(tmp1,tmp2),tmp3);
      theta[eN] = b - tmp4*(a+b);
      theta[eN] = fmax(theta[eN],0.0);

      /*overwrite mass matrix with just correction and element residual*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = (nDOF_trial_element-1)*theta[eN];
	  elementResidual[eN*nDOF_test_element + i] += (nDOF_trial_element-1)*theta[eN]*u_dof[l2g[eN*nDOF_trial_element+i]];

	  for (j=0; j < i; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	  for (j=i+1; j < nDOF_trial_element; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	}  
    }/*eN*/


}
extern "C"
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
					       double* slumpedMassMatrixCorrection)
{
  /*compute slumping in 1d assuming only monotinicity constraints on local element but
    some global coupling through delta R (use dr_i = r_i - r_i-1, for rhs)*/
  int eN,i,j,k,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int I0,Im1,Ip1;
  double drm1,dr0,a,b,tmp1,tmp2,tmp3,tmp4;
  int eN_neighbor0;
  assert(nDOF_test_element == 2);
  for (eN=0; eN<nElements_global; eN++)
    {
      I0 = l2g[eN*nDOF_test_element + 0];
      Ip1= l2g[eN*nDOF_test_element + 1];
      dr0 = fabs(rhs[Ip1]-rhs[I0]);
      
      /*local stencil relative to node I*/
      /*I0 neighbor across from node I1 and vice versa*/
      eN_neighbor0 =elementNeighborsArray[eN*nElementBoundaries_element+1];
      

      
      drm1 = 0.0; Im1 = -1;
      if (eN_neighbor0 >= 0)
	{
	  for (j=0; j < nDOF_test_element; j++)
	    {
	      if (l2g[eN_neighbor0*nDOF_test_element + j] != I0)
		{
		  Im1 = l2g[eN_neighbor0*nDOF_test_element + j];
		  drm1 = fabs(rhs[I0]-rhs[Im1]);
		  break;
		}
	    }
	}
      /*calculate local mass matrix entries isn't exactly right
	for now, since assumes matrix same on both elements*/
      a = 0.0; b = 0.0; i = 0; j = 1;
      for (k=0; k < nQuadraturePoints_element; k++)
	{
	  a += 
	    dm[eN*nQuadraturePoints_element+k]
	    *
	    v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
	      i]
	    *
	    w[eN*nQuadraturePoints_element*nDOF_trial_element + 
	      k*nDOF_test_element + 
	      i]
	    *
	    dV[eN*nQuadraturePoints_element+k];
	  b += 
	    dm[eN*nQuadraturePoints_element+k]
	    *
	    v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
	      j]
	    *
	    w[eN*nQuadraturePoints_element*nDOF_trial_element + 
	      k*nDOF_test_element + 
	      i]
	    *
	    dV[eN*nQuadraturePoints_element+k];
	}
      /*mwf what if account for coupling in mass matrix entries only?
	a *= 2.0;
      */
      tmp1 = fmin(dr0,drm1);
      tmp2 = tmp1/(dr0+drm1+1.0e-12);
      tmp3 = b-(a+b)*tmp2;
      tmp4 = fmax(tmp3,0.0);
      theta[eN] = tmp4;

      /*overwrite mass matrix with just correction and element residual*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = theta[eN];
	  elementResidual[eN*nDOF_test_element + i] += theta[eN]*u_dof[l2g[eN*nDOF_trial_element+i]];

	  for (j=0; j < i; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	  for (j=i+1; j < nDOF_trial_element; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	}  
    }/*eN*/


}

extern "C"
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
					 double* slumpedMassMatrixCorrection)
{
  /*
    try stupid approach where just slump local element mass matrix
   */
  int eN,i,j,k,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int I0,I1,I2;
  double dr0,dr1,dr2,a,b,tmp1,tmp2,tmp3,tmp4;
  
  assert(nDOF_test_element == 3);
  for (eN=0; eN<nElements_global; eN++)
    {
      I0 = l2g[eN*nDOF_test_element + 0];
      I1 = l2g[eN*nDOF_test_element + 1];
      I2 = l2g[eN*nDOF_test_element + 2];
      dr0 = fabs(rhs[I1]-rhs[I0]);
      dr1 = fabs(rhs[I2]-rhs[I1]);
      dr2 = fabs(rhs[I0]-rhs[I2]);
      /*calculate local mass matrix, and store in correction for now */
      for (i=0; i < nDOF_test_element; i++)
	for (j=0; j < nDOF_trial_element; j++)
	  {
	    slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = 0.0;
	    for (k=0; k < nQuadraturePoints_element; k++)
	      {
		slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] +=  
		  dm[eN*nQuadraturePoints_element+k]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
		    j]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+k];
	      }
	  }
      a = 0.0; b=123456.0; /*take max,min to make sure don't have some small roundoff?*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  a = fmax(a,slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i]);
	  for (j=0; j < nDOF_trial_element; j++)
	    {
	      b = fmin(b,slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j]);
	    }
	}

      b = fmax(b,0.0); a = fmax(a,0.0);
      tmp1 = fmin(dr0,dr1); tmp1 = fmin(tmp1,dr2);
      tmp2 = tmp1*(a+b)/(dr0+dr1+dr2+1.0e-12);
      tmp3 = b-tmp2;
      theta[eN] = fmax(tmp3,0.0);
      /*mwf debug*/
      if (tmp1 > 1.0e-5)
	{
	  std::cout<<"slump2d eN="<<eN<<" a="<<a<<" b="<<b<<"dr=["<<dr0<<","<<dr1<<","<<dr2<<"]"<<std::endl;
	  /*mwf debug*/
	  std::cout<<"\t min(dr)="<<tmp1<<" min(dr)*(a+b)/(r0+r1+r2)="<<tmp2<<" b-tmp2="<<tmp3<<" theta="<<theta[eN]<<std::endl;
	}
      /*update mass matrix, storing only correction and element residual*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = (nDOF_trial_element-1)*theta[eN];
	  elementResidual[eN*nDOF_test_element + i] += (nDOF_trial_element-1)*theta[eN]*u_dof[l2g[eN*nDOF_trial_element+i]];

	  for (j=0; j < i; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	  for (j=i+1; j < nDOF_trial_element; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	}  
    }/*eN*/
}
extern "C"
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
					 double* slumpedMassMatrixCorrection)
{
  /*
    try stupid approach where just slump local element mass matrix
   */
  int eN,i,j,k,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int I0,I1,I2;
  double dr0,dr1,dr2,a,b,tmp0,tmp1,tmp2,tmp3,tmp4;
  assert(nDOF_test_element == 3);
  int debugLocalMassMatrix = 0;
  for (eN=0; eN<nElements_global; eN++)
    {
      I0 = l2g[eN*nDOF_test_element + 0];
      I1 = l2g[eN*nDOF_test_element + 1];
      I2 = l2g[eN*nDOF_test_element + 2];

      dr0 = fabs(rhs[I1]-rhs[I0]);
      dr1 = fabs(rhs[I2]-rhs[I1]);
      dr2 = fabs(rhs[I0]-rhs[I2]);

      /*calculate local mass matrix, and store in correction for now */
      for (i=0; i < nDOF_test_element; i++)
	for (j=0; j < nDOF_trial_element; j++)
	  {
	    slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = 0.0;
	    for (k=0; k < nQuadraturePoints_element; k++)
	      {
		slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] +=  
		  dm[eN*nQuadraturePoints_element+k]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
		    j]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+k];
	      }
	  }
      a = 0.0; b=123456.0; /*take max,min to make sure don't have some small roundoff?*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  a = fmax(a,slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i]);
	  for (j=0; j < nDOF_trial_element; j++)
	    {
	      b = fmin(b,slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j]);
	    }
	}

      b = fmax(b,0.0); a = fmax(a,0.0);
      /*arbitrarily scaling matrices to account for coupling in mass matrix entries doesn't seem to make a big difference?
	eg in a structured 2d mesh b *= 2.0; a *= 6.0;
      */
      tmp0 = (b*(dr1+dr2) - (a+b)*dr0)/(dr0+dr1+dr2+1.0e-12);
      tmp1 = (b*(dr0+dr2) - (a+b)*dr1)/(dr0+dr1+dr2+1.0e-12);
      tmp2 = (b*(dr0+dr1) - (a+b)*dr2)/(dr0+dr1+dr2+1.0e-12);
      tmp3 = fmax(tmp0,tmp1); tmp3 = fmax(tmp2,tmp3);
      tmp3 *= adjustFactor;/*add a safety or sharpening factor?*/
      theta[eN] = fmax(tmp3,0.0);
      /*mwf debug*/
      if (tmp1 > 1.0e-5 && false)
	{
	  std::cout<<"slump2d eN="<<eN<<" a="<<a<<" b="<<b<<"dr=["<<dr0<<","<<dr1<<","<<dr2<<"]"<<std::endl;
	  /*mwf debug*/
	  std::cout<<"\t min(dr)="<<tmp1<<" min(dr)*(a+b)/(r0+r1+r2)="<<tmp2<<" b-tmp2="<<tmp3<<" theta="<<theta[eN]<<std::endl;
	}
      /*mwf debug test that I am getting a positive lhs?*/
      if (debugLocalMassMatrix)
	{
	  /*build actual matrix*/
	  for (i=0; i < nDOF_test_element; i++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = a + (nDOF_trial_element-1)*theta[eN];
	      for (j=0; j < i; j++)
		{
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = b-theta[eN];
		}
	      for (j=i+1; j < nDOF_trial_element; j++)
		{
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = b-theta[eN];
		}
	    }
	  double rhs_test[3] = {dr0,dr1,dr2}; double scale=1.0;
	  PROTEUS_LAPACK_INTEGER pivots[3] = {0,0,0}; PROTEUS_LAPACK_INTEGER pivotsCol[3] = {0,0,0};
	  PROTEUS_LAPACK_INTEGER info=0;PROTEUS_LAPACK_INTEGER nsys=3;
	  //std::cout<<"in Slump2d rhs = [";
	  //for (i=0; i < 3; i++)
	  //{
	  //std::cout<<" "<<rhs_test[i];
	  //}
	  //std::cout<<"]"<<std::endl;
	  dgetc2_(&nsys,
		  &slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element],
		  &nsys,
		  pivots,
		  pivotsCol,
		  &info);
	  dgesc2_(&nsys,
		  &slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element],
		  &nsys,
		  rhs_test,
		  pivots,
		  pivotsCol,
		  &scale);
	  //std::cout<<"solution = [";
	  for (i=0; i < 3; i++)
	    {
	      //std::cout<<" "<<rhs_test[i];
	      assert(rhs_test[i] > -1.0e-6);
	    }
	  //std::cout<<"]"<<std::endl;
	  /*cleanup*/
	  for (i=0; i < nDOF_test_element; i++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = 0.0;
	      for (j=0; j < i; j++)
		{
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = 0.0;
		}
	      for (j=i+1; j < nDOF_trial_element; j++)
		{
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = 0.0;
		}
	    }

	}
      /*update mass matrix, storing only correction and element residual*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = (nDOF_trial_element-1)*theta[eN];
	  elementResidual[eN*nDOF_test_element + i] += (nDOF_trial_element-1)*theta[eN]*u_dof[l2g[eN*nDOF_trial_element+i]];

	  for (j=0; j < i; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	  for (j=i+1; j < nDOF_trial_element; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	}  
    }/*eN*/
}

/*attempt to extend Tom's approach to 2d by looking at a 1d tri-diagonal system in the 'upwind' direction?*/
extern "C"
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
						const int* nodeStarOffsets,
						const int* nodeStarArray,
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
						double* slumpedMassMatrixCorrection)
{
  int nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  const int nSpace = 2;
  double volume,df_avg[3] = {0.0,0.0,0.0}, edge_IJ[3] = {0.,0.0,.0,}, edge_Ip1J[3] = {0.,0.,0.};
  assert(nDOF_test_element == 3);
  for (int eN=0; eN<nElements_global; eN++)
    {

      //compute average 'velocity' for upwinding
      for (int J=0; J < nSpace; J++)
	  df_avg[J] = 0.0;
      volume = 0.0;
      for (int k=0; k < nQuadraturePoints_element; k++)
	{
	  volume += dV[eN*nQuadraturePoints_element + k];
	  for (int J=0; J < nSpace; J++)
	    {
	      df_avg[J] += df[eN*nQuadraturePoints_element*nSpace + k*nSpace +J]*
		dV[eN*nQuadraturePoints_element + k];
	    }
	}
      for (int J=0; J < nSpace; J++)
	  df_avg[J] /= volume;
      
      //find the most 'upwind' node and call it I0
      double max_df_ni = 0.0; int i_max_df_ni=0;
      for (int i=0; i < nNodes_element; i++)
	{
	  double df_ni = 0.0;
	  //dot product with outer normal across from i
	  for (int J=0; J < nSpace; J++)
	    df_ni += df_avg[J]*elementBoundaryLocalOuterNormalsArray[eN*nElementBoundaries_element*nSpace + i*nSpace + J];
	  if (df_ni > max_df_ni)
	    {
	      i_max_df_ni = i; max_df_ni = df_ni;
	    } 
	}
      int I0 = l2g[eN*nDOF_test_element + i_max_df_ni];

      //would like to extract a local tridiagonal system centered around I0 that is as aligned as possible with the
      //dominant flow direction

      //loop through nodes {J} in I0 node star and pick out the most 'upwind' and 'downwind' neighboring nodes
      // call the most upwind node Im1, and the most downwind one Ip1
      //To determine alignment with velocity, look at edge (normalized by edge length) dotted with df_avg
      
      double max_delta_IJ = 0.0;
      assert(nodeStarOffsets[I0+1] - nodeStarOffsets[I0] > 1);
      int offset_max_delta_IJ  = nodeStarOffsets[I0]; 
      for (int offset = nodeStarOffsets[I0]; offset < nodeStarOffsets[I0+1]; offset++)
	{
	  const int J = nodeStarArray[offset];
	  assert(J >= 0);
	  assert(J < nNodes_global);

	  double length = 0.0;
	  for (int l = 0; l < nSpace; l++)
	    {
	      edge_IJ[l] = (nodeArray[J*3+l]-nodeArray[I0*3+l]);
	      length += (nodeArray[J*3+l]-nodeArray[I0*3+l])*(nodeArray[J*3+l]-nodeArray[I0*3+l]);
	    }
	  length = sqrt(length);
	  assert(length  > 0.0);
	  for (int l = 0; l < nSpace; l++)
	    edge_IJ[l] /= length;
	  double delta_IJ = 0.0;
	  for (int l = 0; l < nSpace; l++)
	    delta_IJ += edge_IJ[l]*df_avg[l];
	  if (delta_IJ > max_delta_IJ || offset == nodeStarOffsets[I0])
	    {
	      offset_max_delta_IJ = offset;
	      max_delta_IJ = delta_IJ;
	    }
	}

      double min_delta_IJ = 0.0; int offset_min_delta_IJ = nodeStarOffsets[I0+1]-1; 
      for (int offset = nodeStarOffsets[I0+1]-1; offset >= nodeStarOffsets[I0]; offset--)
	{
	  const int J = nodeStarArray[offset];
	  assert(J >= 0);
	  assert(J < nNodes_global);

	  double length = 0.0;
	  for (int l = 0; l < nSpace; l++)
	    {
	      edge_IJ[l] = (nodeArray[J*3+l]-nodeArray[I0*3+l]);
	      length += (nodeArray[J*3+l]-nodeArray[I0*3+l])*(nodeArray[J*3+l]-nodeArray[I0*3+l]);
	    }
	  length = sqrt(length);
	  assert(length  > 0.0);
	  for (int l = 0; l < nSpace; l++)
	    edge_IJ[l] /= length;
	  double delta_IJ = 0.0;
	  for (int l = 0; l < nSpace; l++)
	    delta_IJ += edge_IJ[l]*df_avg[l];
	  if (delta_IJ < min_delta_IJ || offset == nodeStarOffsets[I0+1]-1)
	    {
	      offset_min_delta_IJ = offset;
	      min_delta_IJ = delta_IJ;
	    }
	}
      //neighbor with most 'aligned' edge is downwind since edge defined as J-I0
      //at a corner, say might get right angled configuration here if just take raw min/max
      int Ip1= nodeStarArray[offset_max_delta_IJ];
      int Im1= nodeStarArray[offset_min_delta_IJ];
      assert(Ip1 != Im1);
      double dr0  = fabs(rhs[Ip1]-rhs[I0]);
      double drm1 = fabs(rhs[I0]-rhs[Im1]);
 
      //now need to find most downwind neighbor in Ip1's node star (down-down wind)
      double max_delta_Ip1J = 0.0;
      assert(nodeStarOffsets[Ip1+1] - nodeStarOffsets[Ip1]> 1);
      int offset_max_delta_Ip1J = nodeStarOffsets[Ip1];
      for (int offset = nodeStarOffsets[Ip1]; offset < nodeStarOffsets[Ip1+1]; offset++)
	{
	  const int J = nodeStarArray[offset];
	  assert(J >= 0);
	  assert(J < nNodes_global);
	  
	  double length = 0.0;
	  for (int l = 0; l < nSpace; l++)
	    {
	      edge_Ip1J[l] = (nodeArray[J*3+l]-nodeArray[Ip1*3+l]);
	      length += (nodeArray[J*3+l]-nodeArray[Ip1*3+l])*(nodeArray[J*3+l]-nodeArray[Ip1*3+l]);
	    }
	  length = sqrt(length);
	  assert(length  > 0.0);
	  for (int l = 0; l < nSpace; l++)
	    edge_Ip1J[l] /= length;
	  double delta_Ip1J = 0.0;
	  for (int l = 0; l < nSpace; l++)
	    delta_Ip1J += edge_Ip1J[l]*df_avg[l];
	  if (delta_Ip1J > max_delta_Ip1J || offset == nodeStarOffsets[Ip1])
	    {
	      offset_max_delta_Ip1J = offset;
	      max_delta_Ip1J = delta_Ip1J;
	    }
	}
      int Ip2 = -1; double drp1 = 0.0; 
      if (nodeStarArray[offset_max_delta_Ip1J] != I0 && 
	  nodeStarArray[offset_max_delta_Ip1J] != Im1)
	Ip2 = nodeStarArray[offset_max_delta_Ip1J];
	  	  
     
      if (Ip2 >= 0)
	drp1 = fabs(rhs[Ip2]-rhs[Ip1]);
      
      /*calculate local mass matrix entries isn't exactly right
	for now, since assumes matrix same on all elements*/
      double a = 0.0; double b = 0.0; int i = 0; int j = 1;
      for (int k=0; k < nQuadraturePoints_element; k++)
	{
	  a += 
	    dm[eN*nQuadraturePoints_element+k]
	    *
	    v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
	      i]
	    *
	    w[eN*nQuadraturePoints_element*nDOF_trial_element + 
	      k*nDOF_test_element + 
	      i]
	    *
	    dV[eN*nQuadraturePoints_element+k];
	  b += 
	    dm[eN*nQuadraturePoints_element+k]
	    *
	    v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
	      j]
	    *
	    w[eN*nQuadraturePoints_element*nDOF_trial_element + 
	      k*nDOF_test_element + 
	      i]
	    *
	    dV[eN*nQuadraturePoints_element+k];
	}
      
      double tmp1 = dr0/(drm1+dr0+drp1+1.0e-12);
      double tmp2 = fmin(drm1,dr0)/(drm1+dr0+1.0e-12);
      double tmp3 = fmin(drp1,dr0)/(drp1+dr0+1.0e-12);
      double tmp4 = fmin(fmin(tmp1,tmp2),tmp3);
      theta[eN] = b - tmp4*(a+b);
      theta[eN] = fmax(theta[eN],0.0);

      /*overwrite mass matrix with just correction and element residual*/
      for (int i=0; i < nDOF_test_element; i++)
	{
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = (nDOF_trial_element-1)*theta[eN];
	  elementResidual[eN*nDOF_test_element + i] += (nDOF_trial_element-1)*theta[eN]*u_dof[l2g[eN*nDOF_trial_element+i]];

	  for (int j=0; j < i; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	  for (int j=i+1; j < nDOF_trial_element; j++)
	    {
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	}  
    }/*eN*/


}

extern "C"
void updateElementJacobianWithSlumpedMassApproximation(int nElements_global,
						       int nDOF_trial_element,
						       int nDOF_test_element,
						       const double* theta,
						       double* elementJacobianMatrix)
{
  int eN,i,j,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  for (eN = 0; eN < nElements_global; eN++)
    {
      for (i=0; i < nDOF_test_element; i++)
	{
	  for (j=0; j < i; j++)
	    elementJacobianMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] -= theta[eN];
	  elementJacobianMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] += (nDOF_trial_element-1)*theta[eN];
	  for (j=i+1; j < nDOF_trial_element; j++)
	    elementJacobianMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] -= theta[eN];
	}
    }

}


extern "C"
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
					       double* slumpedMassMatrixCorrection)
{
  int eN,i,j,k,kk,I,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int J,K,eN_neighbor;
  double mii,mij,mik,sei,theij;
  assert(nDOF_test_element == 2);
  for (eN=0; eN<nElements_global; eN++)
    {
      for (i=0; i < nDOF_test_element; i++)
	{
	  j = (i+1) % nDOF_test_element;
	  I = l2g[eN*nDOF_test_element + i];
	  J = l2g[eN*nDOF_test_element + j];
	  /*compute local element mass matrix contributions*/
	  mii = 0.0; mij=0.0;
	  for (kk=0; kk < nQuadraturePoints_element; kk++)
	      {
		mii += 
		  dm[eN*nQuadraturePoints_element+kk]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_trial_element + 
		    i]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+kk];
		mij += 
		  dm[eN*nQuadraturePoints_element+kk]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_trial_element + 
		    j]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+kk];
	      }

	  eN_neighbor = elementNeighborsArray[eN*nElementBoundaries_element+j]; /*across from j*/
      
	  K=-1;theij=0.0; /*lump by default*/
	  if (eN_neighbor >= 0)
	    {
	      for (k=0; k < nDOF_test_element; k++)
		{
		  if (l2g[eN_neighbor*nDOF_test_element + k] != I)
		    {
		      K = l2g[eN_neighbor*nDOF_test_element + k];
		      break;
		    }
		}
	      assert(K >= 0); mik = 0.0;
	      for (kk=0; kk < nQuadraturePoints_element; kk++)
		{
		mik += 
		  dm[eN_neighbor*nQuadraturePoints_element+kk]
		  *
		  v[eN_neighbor*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_trial_element + 
		    k]
		  *
		  w[eN_neighbor*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_test_element + 
		    i]
		  *
		  dV[eN_neighbor*nQuadraturePoints_element+kk];
		}
	      sei = 0.0;
	      if (fabs(u_dof_limit[I]-u_dof_limit[K]) > 1.0e-6)
		sei = (u_dof_limit[J]-u_dof_limit[I])/(u_dof_limit[I]-u_dof_limit[K]);
	      if (sei > 0.0 && fabs(sei) > 1.0e-6)
		{
		  theij = mik/sei - mij;
		  theij = fmax(theij,0.0);
		}
	    }
	  /*mwf hack force lumping
	    theij = 0.0;
	  */
	  /*update residual and modified mass matrix*/
	  elementResidual[eN*nDOF_test_element + i] += (mij + theij)*u_dof[I] - (mij+theij)*u_dof[J];
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = mij + theij;
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -(mij + theij);
	  /*mwf hack try lumped mass matrix as approximate jacobian?
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] =  mij;
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -mij;
	  */
	}/*i loop*/
    }/*eN*/

}
double modifiedVanLeer(double r)
{
  return (fabs(r)+r)/(1. + fmax(1.,fabs(r)));
}

extern "C"
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
						  double* slumpedMassMatrixCorrection)
{
  int eN,i,j,k,kk,I,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int J,K,eN_neighbor;
  double mii,mij,mik,sei,theij,theijmax;
  assert(nDOF_test_element == 2);
  for (eN=0; eN<nElements_global; eN++)
    {
      /*compute only one lumping parameter per element */
      for (i=0; i < nDOF_test_element; i++)
	{
	  j = (i+1) % nDOF_test_element;
	  I = l2g[eN*nDOF_test_element + i];
	  J = l2g[eN*nDOF_test_element + j];
	  /*compute local element mass matrix contributions*/
	  mii = 0.0; mij=0.0;
	  for (kk=0; kk < nQuadraturePoints_element; kk++)
	      {
		mii += 
		  dm[eN*nQuadraturePoints_element+kk]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_trial_element + 
		    i]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+kk];
		mij += 
		  dm[eN*nQuadraturePoints_element+kk]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_trial_element + 
		    j]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+kk];
	      }

	  eN_neighbor = elementNeighborsArray[eN*nElementBoundaries_element+j]; /*across from j*/
      
	  K=-1;theij=0.0; /*lump by default*/
	  if (eN_neighbor >= 0)
	    {
	      for (k=0; k < nDOF_test_element; k++)
		{
		  if (l2g[eN_neighbor*nDOF_test_element + k] != I)
		    {
		      K = l2g[eN_neighbor*nDOF_test_element + k];
		      break;
		    }
		}
	      assert(K >= 0); mik = 0.0;
	      for (kk=0; kk < nQuadraturePoints_element; kk++)
		{
		mik += 
		  dm[eN_neighbor*nQuadraturePoints_element+kk]
		  *
		  v[eN_neighbor*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_trial_element + 
		    k]
		  *
		  w[eN_neighbor*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_test_element + 
		    i]
		  *
		  dV[eN_neighbor*nQuadraturePoints_element+kk];
		}
	      sei = 0.0;
	      if (fabs(u_dof_limit[I]-u_dof_limit[K]) > 1.0e-6)
		sei = (u_dof_limit[J]-u_dof_limit[I])/(u_dof_limit[I]-u_dof_limit[K]);
	      if (sei > 0.0 && fabs(sei) > 1.0e-6)
		{
		  theij = mik/sei - mij;
		  theij = fmax(theij,0.0);/*modifiedVanLeer(theij);fmax(theij,0.0);*/
		}
	    }
	  //theijmax += theij/nDOF_test_element; /*fmax(theijmax,theij);*/
	  if (i == 0) 
	    theijmax = theij;
	  else
	    theijmax = fmin(theij,theijmax);
	  //try to put a lower bound on the diffusion?
	  theijmax = fmin(theijmax,mij);
	}/*first loop to pick theta*/
      i = 0; j = 1;
      I = l2g[eN*nDOF_test_element + i];
      J = l2g[eN*nDOF_test_element + j];
      
      /*update residual and modified mass matrix*/
      elementResidual[eN*nDOF_test_element + i] += (mij + theijmax)*u_dof[I] - (mij+theijmax)*u_dof[J];
      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = mij + theijmax;
      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -(mij + theijmax);
      elementResidual[eN*nDOF_test_element + j] += (mij + theijmax)*u_dof[J] - (mij+theijmax)*u_dof[I];
      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + j*nDOF_trial_element + j] = mij + theijmax;
      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + j*nDOF_trial_element + i] = -(mij + theijmax);
      /*mwf hack try lumped mass matrix as approximate jacobian?
	slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] =  mij;
	slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -mij;
      */
    }/*eN*/

}


extern "C"
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
					       double* slumpedMassMatrixCorrection)
{
  int eN,i,j,k,kk,I,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int J,K,eN_neighbor;
  double mii,mij,mik,sei,rei,theij,theik;
  assert(nDOF_test_element == 3);
  for (eN=0; eN<nElements_global; eN++)
    {
      for (i=0; i < nDOF_test_element; i++)
	{
	  j = (i+1) % nDOF_test_element;
	  k = (i+2) % nDOF_test_element;
	  I = l2g[eN*nDOF_test_element + i];
	  J = l2g[eN*nDOF_test_element + j];
	  K = l2g[eN*nDOF_test_element + k];
	  /*compute local element mass matrix contributions*/
	  mii = 0.0; mij=0.0; mik=0.0;
	  for (kk=0; kk < nQuadraturePoints_element; kk++)
	    {
	      mii += 
		  dm[eN*nQuadraturePoints_element+kk]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_trial_element + 
		    i]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    kk*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+kk];
	      mij += 
		dm[eN*nQuadraturePoints_element+kk]
		*
		v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		  kk*nDOF_trial_element + 
		  j]
		*
		w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		  kk*nDOF_test_element + 
		  i]
		*
		dV[eN*nQuadraturePoints_element+kk];
	      mik += 
		dm[eN*nQuadraturePoints_element+kk]
		*
		v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		  kk*nDOF_trial_element + 
		  k]
		*
		w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		  kk*nDOF_test_element + 
		  i]
		*
		dV[eN*nQuadraturePoints_element+kk];
	    }

	      
	  theij=0.0; theik=0.0;/*lump by default*/

	  sei = 0.0; rei = 0.0;
	  if (fabs(u_dof_limit[J]-u_dof_limit[I]) > 1.0e-6)
	    sei = (u_dof_limit[I]-u_dof_limit[K])/(u_dof_limit[J]-u_dof_limit[I]);
	  if (fabs(sei) > 1.0e-6)
	    rei = 1.0/sei;
	  if (sei > 0.0)
	    {
	      theij = mik*sei - mij;
	      theij = fmax(theij,0.0);
	    }
	  if (rei > 0.0)
	    {
	      theik = mij*rei - mik;
	      theik = fmax(theik,0.0);
	    }
	    
	  /*mwf hack force lumping
	  theij = 0.0; theik = 0.0;
	  */

	  /*update residual and modified mass matrix*/
	  elementResidual[eN*nDOF_test_element + i] += (mij + mik + theij + theik)*u_dof[I] - (mij+theij)*u_dof[J] - (mik + theik)*u_dof[K];
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = mij + theij + mik + theik;
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -(mij + theij);
	  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + k] = -(mik + theik);
	}/*i loop*/
    }/*eN*/

}


extern "C"
void updateElementJacobianWithSlumpedMassCorrection(int nElements_global,
						    int nDOF_trial_element,
						    int nDOF_test_element,
						    const double* elementMassMatrixCorrection,
						    double* elementJacobianMatrix)
{
  int eN,i,j,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  for (eN = 0; eN < nElements_global; eN++)
    {
      for (i=0; i < nDOF_test_element; i++)
	for (j=0; j < nDOF_test_element; j++)
	  elementJacobianMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] += 
	    elementMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j];
    }
}


//assumes globalMassMatrix has already been zeroed
extern "C"
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
				    double* globalMassMatrix)
{
  for (int eN=0; eN<nElements_global; eN++)
    {
      for (int i=0; i < nDOF_test_element; i++)
	{
	  const int I = l2g[eN*nDOF_test_element + i];
	  for (int j=0; j < nDOF_trial_element; j++)
	    {
	      const int J = l2g[eN*nDOF_trial_element + j];
	      double M_ij = 0.0;
	      for (int k=0; k < nQuadraturePoints_element; k++)
		{
		  M_ij += 
		    dm[eN*nQuadraturePoints_element+k]
		    *
		    v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		      k*nDOF_trial_element + 
		      j]
		    *
		    w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		      k*nDOF_test_element + 
		      i]
		    *
		    dV[eN*nQuadraturePoints_element+k];
		}
	      //can't assume colind is in assending order?
	      for (int m=rowptr[I]; m < rowptr[I+1]; m++)
		{
		  if (J == colind[m])
		    {
		      globalMassMatrix[m] += M_ij;
		      break;
		    }
		}
	    }//j
	}//i
    }//eN
}


extern "C"
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
								   double* slumpedMassMatrixCorrection)
{

  int nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int * alpha_ij = new int[nDOF_test_X_trial_element];
  int forceSymmetry=1;

  for (int eN=0; eN<nElements_global; eN++)
    {
      /*calculate local mass matrix, and store in correction for now */
      for (int i=0; i < nDOF_test_element; i++)
	for (int j=0; j < nDOF_trial_element; j++)
	  {
	    slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = 0.0;
	    for (int k=0; k < nQuadraturePoints_element; k++)
	      {
		slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] +=  
		  dm[eN*nQuadraturePoints_element+k]
		  *
		  v[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_trial_element + 
		    j]
		  *
		  w[eN*nQuadraturePoints_element*nDOF_trial_element + 
		    k*nDOF_test_element + 
		    i]
		  *
		  dV[eN*nQuadraturePoints_element+k];
	      }
	  }
      /*while debugging generate a simple local lookup for i,j to offsets in edgeSlumpingParameters*/
      for (int i=0; i < nDOF_trial_element; i++)
	{
	  const int I = l2g[eN*nDOF_trial_element + i];
	  for (int j=0; j < nDOF_trial_element; j++) //can't assume colind is in assending order?
	    {
	      const int J = l2g[eN*nDOF_trial_element + j];
	      for (int m=rowptr[I]; m < rowptr[I+1]; m++)
		{
		  if (J == colind[m])
		    {
		      alpha_ij[i*nDOF_trial_element + j] = globalEdgeSlumpingParameter[m];
		      break;
		    }
		}
	    }
	}
      /*update mass matrix, storing only correction and element residual*/
      if (forceSymmetry)
	{
	  for (int i=0; i < nDOF_test_element; i++)
	    {
	      //diagonal is M_{ii} = M^L_{i} - \sum_{j\ne i} alpha_ij M^c_{ij}
	      //                   = M^c_{ii} + \sum_{j\ne i} M^c_{ij} - \sum_{j\ne i} alpha_ij M^c_{ij}
	      //so correction is 
	      //      \delta M_{ii}= \sum_{j\ne i} (1- alpha_ij) M^c_{ij}
	  
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = 0.0;
	      //upper half
	      for (int j=i+1; j < nDOF_trial_element; j++)
		{
		  const double theta = (1.0-alpha_ij[i*nDOF_trial_element+j])*slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] += theta;
		  elementResidual[eN*nDOF_test_element + i] += theta*u_dof[l2g[eN*nDOF_trial_element+i]];
		}
	      //offdiagonals are 
	      //            M_{ij} = alpha_{ij} M^c_{ij}
	      //so correction is
	      //     \delta M_{ij} = -M^c_{ij} + alpha_{ij}M^c_{ij}
	      //                   = -(1 - alpha_{ij})M^c_{ij}
	      for (int j=i+1; j < nDOF_trial_element; j++)
		{
		  const double theta = (1.0-alpha_ij[i*nDOF_trial_element+j])*slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta;
		  elementResidual[eN*nDOF_test_element + i] -= theta*u_dof[l2g[eN*nDOF_trial_element+j]];
		}
	    }//end i upper half
	  for (int i=0; i < nDOF_test_element; i++)
	    {
	      //diagonal is M_{ii} = M^L_{i} - \sum_{j\ne i} alpha_ij M^c_{ij}
	      //                   = M^c_{ii} + \sum_{j\ne i} M^c_{ij} - \sum_{j\ne i} alpha_ij M^c_{ij}
	      //so correction is 
	      //      \delta M_{ii}= \sum_{j\ne i} (1- alpha_ij) M^c_{ij}
	  
	      
	      //lower half
	      for (int j=0; j < i; j++)
		{
		  const double theta = -slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + j*nDOF_trial_element + i];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] += theta;
		  elementResidual[eN*nDOF_test_element + i] += theta*u_dof[l2g[eN*nDOF_trial_element+i]];
		}
	      //offdiagonals are 
	      //            M_{ij} = alpha_{ij} M^c_{ij}
	      //so correction is
	      //     \delta M_{ij} = -M^c_{ij} + alpha_{ij}M^c_{ij}
	      //                   = -(1 - alpha_{ij})M^c_{ij}
	      for (int j=0; j < i; j++)
		{
		  const double theta = -slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + j*nDOF_trial_element + i];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta;
		  elementResidual[eN*nDOF_test_element + i] -= theta*u_dof[l2g[eN*nDOF_trial_element+j]];
		}
	    }//end i upper half
	}//end force Symmetry
      else
	{
	  for (int i=0; i < nDOF_test_element; i++)
	    {
	      //diagonal is M_{ii} = M^L_{i} - \sum_{j\ne i} alpha_ij M^c_{ij}
	      //                   = M^c_{ii} + \sum_{j\ne i} M^c_{ij} - \sum_{j\ne i} alpha_ij M^c_{ij}
	      //so correction is 
	      //      \delta M_{ii}= \sum_{j\ne i} (1- alpha_ij) M^c_{ij}
	  
	      slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] = 0.0;
	      for (int j=0; j < i; j++)
		{
		  //equivalent of theta in Tom's approach
		  const double theta = (1.0-alpha_ij[i*nDOF_trial_element+j])*slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] += theta;
		  elementResidual[eN*nDOF_test_element + i] += theta*u_dof[l2g[eN*nDOF_trial_element+i]];
		}
	      for (int j=i+1; j < nDOF_trial_element; j++)
		{
		  const double theta = (1.0-alpha_ij[i*nDOF_trial_element+j])*slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] += theta;
		  elementResidual[eN*nDOF_test_element + i] += theta*u_dof[l2g[eN*nDOF_trial_element+i]];
		}
	      //offdiagonals are 
	      //            M_{ij} = alpha_{ij} M^c_{ij}
	      //so correction is
	      //     \delta M_{ij} = -M^c_{ij} + alpha_{ij}M^c_{ij}
	      //                   = -(1 - alpha_{ij})M^c_{ij}
	      for (int j=0; j < i; j++)
		{
		  const double theta = (1.0-alpha_ij[i*nDOF_trial_element+j])*slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta;
		  elementResidual[eN*nDOF_test_element + i] -= theta*u_dof[l2g[eN*nDOF_trial_element+j]];
		}
	      for (int j=i+1; j < nDOF_trial_element; j++)
		{
		  const double theta = (1.0-alpha_ij[i*nDOF_trial_element+j])*slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j];
		  slumpedMassMatrixCorrection[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = -theta;
		  elementResidual[eN*nDOF_test_element + i] -= theta*u_dof[l2g[eN*nDOF_trial_element+j]];
		}
	    }  

	}

    }/*eN*/
  //cleanup
  delete [] alpha_ij;
}


extern "C"
void computeSlumpingParametersFCT_KuzminMoeller10(const int nDOF_global,
						  const int* rowptr, //sparsity information
						  const int* colind,
						  const double* u_dof,
						  const double* u_dof_limit,
						  const double* Mc, //assembled global mass matrix
						  double * Rip,  //waste some temporary space until get 
						  double * Rim,  //algorithm debugged
						  double* edgeSlumpingParameter) //spare rep for slumping factor for edge ij (i.e. alpha_ij)
{
  /***********************************************************************
  try to apply FCT limiting approach for mass matrix
    M^c * u = M^l * u + (M^c - M^l) * u 
  or writing as a flux, for equation i
    (M^c * u)_i = M^l_i u_i + \sum_{j\ne i} M^c_{ij}(u_j - u_i)
  or 
    (M^c * u)_i = M^l_i u_i + \sum_{j\ne i} F_{ij}
    F_{ij} = M^c_{ij}(u_j - u_i), F_{ji} = - F_{ij}
  
    
  Want to replace F_{ij} by \alpha_{ij}(u^*)F^{ij}
  
  Where u^* is some approximation to the solution (either from the last 
    iterate, a mass-lumped low-order solution, or a backtracked approximation

  First we'll try to use the FCT approach from Kuzmin Moeller etal 2010 JCP
    to compute alpha_ij. Note alpha_ij = alpha_ji \in [0,1] always

  Limiting formulas written for F_{ij} on rhs of equation, so use negative
    of formulas above
  ***********************************************************************/
  //todo zero corrections?
  int forceAbsoluteMaxMin=0;
  int debugSymmetry=0;
  int debugOutput=0;
  //mwf debug test max-min
  const double u_max_absolute = 1.0; const double u_min_absolute = 0.0;
  for (int I=0; I < nDOF_global; I++)
    {
      //determine max and min for node star, for now assume sparsity 
      //gives you what you need (no need to get nodeStarOffsets etc)
      double u_max(0.),u_min(0.),u_I(0.0),u_J(0.);
      u_I =u_dof_limit[I];
      if (forceAbsoluteMaxMin)
	{
	  u_I = fmin(fmax(u_I,u_min_absolute),u_max_absolute);
	}
      for (int m=rowptr[I]; m < rowptr[I+1]; m++)
	{
	  const int J = colind[m];
	  u_J =u_dof_limit[J];
	  if (forceAbsoluteMaxMin)
	    {
	      u_J = fmin(fmax(u_J,u_min_absolute),u_max_absolute);
	    }
	  
	  if (m==rowptr[I] || u_J > u_max)
	    u_max =u_J;
	  if (m==rowptr[I] || u_J < u_min)
	    u_min = u_J;
	}
      //distance to local extremum
      double Qip=u_max - u_I;
      double Qim=u_min - u_I;

      //sum of positive, negative anti-diffusive fluxes
      double Pip(0.),Pim(0.);
      //lumped mass coefficient
      double MiL = 0.0;
      for (int m=rowptr[I]; m < rowptr[I+1]; m++)
	{
	  const int J = colind[m];
	  MiL += Mc[m];
	  if (J != I)
	    {
	      u_J =u_dof_limit[J];
	      if (forceAbsoluteMaxMin)
		{
		  u_J = fmin(fmax(u_J,u_min_absolute),u_max_absolute);
		}
	      
	      double Fij=-Mc[m]*(u_dof_limit[J]-u_I);
	      Pip += fmax(0.0,Fij); 
	      Pim += fmin(0.0,Fij);
	    }
	}

      //Nodal correction factors
      //what should be default if Pi=0, Pim=0?
      double rrip(1.0),rrim(1.0);
      //double rrip(0.0),rrim(0.0);
      if (Pip > 1.0e-15)
	rrip = fmin(1.0,MiL*Qip/Pip);
      if (Pim < -1.0e-15)
	rrim = fmin(1.0,MiL*Qim/Pim);

      //waste space until done debugging 
      Rip[I] = rrip; Rim[I] = rrim;
      if (debugOutput > 2 && (fabs(rrip) < 1.0 || fabs(rrim) < 1.0))
	{
	  std::cout<<"KuzminMoeller I= "<<I<<" Rip= "<<rrip<<" Rim= "<<rrim
		   <<" u_max= "<<u_max<<" u_min= "<<u_min<<" u_I= "<<u_dof_limit[I]
		   <<" MiL= "<<MiL<<" Qip= "<<Qip<<" Qim= "<<Qim<<" Pip= "<<Pip
		   <<" Pim= "<<Pim<<std::endl;

	}
    }// First I loop to set Rip, Rim
  //loop through again and set alpha_ij and enforce symmetry
  for (int I=0; I < nDOF_global; I++)
    {
      double u_I(0.0),u_J(0.);
      //mwf stopped here
      u_I =u_dof_limit[I];
      if (forceAbsoluteMaxMin)
	{
	  u_I = fmin(fmax(u_I,u_min_absolute),u_max_absolute);
	}
      //define and enforce alpha_ij=alpha_ji 
      for (int m=rowptr[I]; m < rowptr[I+1]; m++)
	{
	  const int J = colind[m];
	  //if don't explicitly force symmetry
	  // if (J != I) //could set just upper half but then have to loop through offsets for J
	  //   {
	  //     double Fij=-Mc[m]*(u_dof_limit[J]-u_dof_limit[I]);
	  //     
	  //     double Rij = Rip[I];
	  //     if (Fij < 0.0)
	  // 	Rij = Rim[I];
	  //     //j --> i
	  //     double Fji = - Fij;
	  //     double Rji = Rip[J];
	  //     if (Fji < 0.0)
	  // 	Rji = Rim[J];
	  //     double aij = fmin(Rij,Rji);
	  //     //mwf hack force lumping
	  //     //aij = 0.0;
	  //     //test constant
	  //     //aij = 0.5;
	  //     assert(0.0 <= aij);
	  //     assert(1.0 >= aij);
	  //     edgeSlumpingParameter[m] = aij;

	  //   }
	  if (J != I) //could set just upper half but then have to loop through offsets for J
	    {
	      u_J =u_dof_limit[J];
	      if (forceAbsoluteMaxMin)
		{
		  u_J = fmin(fmax(u_J,u_min_absolute),u_max_absolute);
		}

	      double Fij=-Mc[m]*(u_J-u_I);
	      
	      double Rij = Rip[I];
	      if (Fij < 0.0)
		Rij = Rim[I];
	      //j --> i
	      double Fji = -Fij;
	      double Rji = Rip[J];
	      if (Fji < 0.0)
		Rji = Rim[J];
	      double aij = fmin(Rij,Rji);
	      //mwf hack force lumping
	      //aij = 0.0;
	      //test constant
	      //aij = 0.5;
	      assert(0.0 <= aij);
	      assert(1.0 >= aij);
	      edgeSlumpingParameter[m] = aij;
	      for (int mm=rowptr[J]; mm < rowptr[J+1]; mm++)
		{
		  const int II = colind[mm];
		  if (II == I)
		    {
		      edgeSlumpingParameter[mm] = aij;
		      break;
		    }
		}
	    }
	  else
	    edgeSlumpingParameter[m] = 0.0; //what's a better default? 1.0
	}
    }//loop through I to set aj

  //test symmetry
  if (debugSymmetry)
    {
      for (int I=0; I < nDOF_global; I++)
	{
	  for (int m=rowptr[I]; m < rowptr[I+1]; m++)
	    {
	      const int J = colind[m];
	      const double mij = Mc[m];
	      for (int mm = rowptr[J]; mm < rowptr[J+1]; mm++)
		{
		  const int II = colind[mm];
		  if (II == I)
		    {
		      const double mji = Mc[mm];
		      if (fabs(mij-mji) > 1.0e-12)
			{
			  std::cout<<"problem in Kuzmin_Turek limiting I= "<<I<<" J= "<<J<<" mij= "<<mij
				   <<" mji = "<<mji<<std::endl;
			  assert(0);
			}
		    }
		}
	    }
	}
      for (int I=0; I < nDOF_global; I++)
	{
	  for (int m=rowptr[I]; m < rowptr[I+1]; m++)
	    {
	      const int J = colind[m];
	      const double aij = edgeSlumpingParameter[m];
	      for (int mm = rowptr[J]; mm < rowptr[J+1]; mm++)
		{
		  const int II = colind[mm];
		  if (II == I)
		    {
		      const double aji = edgeSlumpingParameter[mm];
		      if (fabs(aij-aji) > 1.0e-12)
			{
			  std::cout<<"problem in Kuzmin_Turek limiting I= "<<I<<" J= "<<J<<" aij= "<<aij
				   <<" aji = "<<aji<<std::endl;
			  assert(0);
			}
		    }
		}
	    }
	}
    }//end debug symmetry
}

extern "C"
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
			    std::vector<double>& x_gq,
			    std::vector<double>& dV_gq)
{
  std::map<int,std::list<double> > elementsToTrackedPoints;

  //determine which elements contain SSIPs
  ELLAM::collectTrackedPointsInElements(nElements_global,
				       nNodes_element,
				       nElementBoundaries_element,
				       nSpace,
				       nPoints_tracked,
				       boundaryTolerance,
				       neighborTolerance,
				       nodeArray,
				       elementNodesArray,
				       elementBoundariesArray,
				       elementBoundaryLocalOuterNormalsArray,
				       elementBoundaryBarycentersArray,
				       element_track,
				       flag_track,
				       x_track,
				       elementsToTrackedPoints);

  //reset these 
  elements_gq.clear(); x_gq.clear(); dV_gq.clear();
  //conservative enough guess on memory needed?
  elements_gq.reserve(nElements_global*nPointsPerElement_default*2);
  dV_gq.reserve(nElements_global*nPointsPerElement_default*2);
  x_gq.reserve(nElements_global*nPointsPerElement_default*2*3);
  
  //loop through elements, 
  //if element has SSIPs, then triangulate element and
  // create new quadrature points on each sub-element
  //if now SSIPs in element, just copy over default element quadrature points
  if (nSpace == 1)
    {
      std::vector<double> x_element,w_element;
      for (int eN=0; eN < nElements_global; eN++)
	{
	  int nPoints_eN = 0;
	  if (elementsToTrackedPoints.find(eN) != elementsToTrackedPoints.end())
	    {
	      ELLAM::createElementSSIPs_1d(elementsToTrackedPoints[eN],
					   x_element,
					   w_element);
	      nPoints_eN = int(w_element.size());
	      for (int k=0; k < nPoints_eN; k++)
		{
		  dV_gq.push_back(w_element[k]);
		  for (int I=0; I < 3; I++)
		    x_gq.push_back(x_element[k*3+I]);
		}
	    }
	  else
	    {
	      nPoints_eN = nPointsPerElement_default;
	      for (int k=0; k < nPointsPerElement_default; k++)
		{
		  const double wk = q_dV_default[eN*nPointsPerElement_default +k];
		  dV_gq.push_back(wk);
		  for (int I=0; I < 3; I++)
		    {
		      const double xkI=q_x_default[eN*nPointsPerElement_default*3 +k*3 +I];
		      x_gq.push_back(xkI);
		    }
		}
	    }//generating points for eN
	  assert(nPoints_eN > 0);
	  for (int k=0; k < nPoints_eN; k++)
	    {
	      elements_gq.push_back(eN);
	    }
	}//eN
    }//1d
  else
    {
      std::cout<<"generateQuadratureArraysForSSIPs nSpace= "<<nSpace<<" not implemented! "<<std::endl;
      assert(0);
    }

}

extern "C"
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
				   std::vector<double>& x_ssip)

{
  std::map<int,std::list<double> > elementsToTrackedPoints;
  assert(element_offsets_ssip.size() == unsigned(nElements_global+1));
  //determine which elements contain SSIPs
  ELLAM::collectTrackedPointsInElements(nElements_global,
				       nNodes_element,
				       nElementBoundaries_element,
				       nSpace,
				       nPoints_tracked,
				       boundaryTolerance,
				       neighborTolerance,
				       nodeArray,
				       elementNodesArray,
				       elementBoundariesArray,
				       elementBoundaryLocalOuterNormalsArray,
				       elementBoundaryBarycentersArray,
				       element_track,
				       flag_track,
				       x_track,
				       elementsToTrackedPoints);

  //reset these 
  x_ssip.clear(); 
  //conservative enough guess on memory needed?
  x_ssip.reserve(nPoints_tracked*3);
  
  //loop through elements, 
  element_offsets_ssip[0] = 0;

  std::vector<double> x_element;
  for (int eN=0; eN < nElements_global; eN++)
    {
      int nPoints_eN = 0;
      if (elementsToTrackedPoints.find(eN) != elementsToTrackedPoints.end())
	{
	  nPoints_eN = int(elementsToTrackedPoints[eN].size()/3);
	  std::list<double>::iterator eN_it = elementsToTrackedPoints[eN].begin();
	  while (eN_it != elementsToTrackedPoints[eN].end())
	    {
	      x_ssip.push_back(*eN_it); 
	      eN_it++;
	    }
	}
      element_offsets_ssip[eN+1] = element_offsets_ssip[eN]+nPoints_eN;
    }

}


/**********************************************************************
 utility routines for SSIPs
 **********************************************************************/
namespace ELLAM
{
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
				   std::map<int,std::list<double> >& elementsToTrackedPoints)
{

  /***********************************************************************
     loop through points, and assign them to element containing them
     if the point is within boundaryTolerance of boundary, project to
     the boundary first. Then add to the element if the point is not 
     too close to another point (using neighborTolerance)
     in the element, including vertices of the mesh

   ***********************************************************************/
  double x[3] = {0.0,0.0,0.0};
  //
  elementsToTrackedPoints.clear();

  for (int k=0; k < nPoints_tracked; k++)
    {
      //copy over for some local manipulations
      x[0] = x_track[k*3+0]; x[1] = x_track[k*3+1]; x[2] = x_track[k*3+2];
      const int eN = element_track[k];
      if (eN >= 0 && flag_track[k] >= -1)
	{
	  //is point too close to the boundary?
	  for (int ebN = 0; ebN < nElementBoundaries_element; ebN++)
	    {
	      const int ebN_global = elementBoundariesArray[eN*nElementBoundaries_element+ebN];
	      double dxf = 0.0;
	      for (int I = 0; I < nSpace; I++)
		{
		  dxf += (elementBoundaryBarycentersArray[ebN_global*3+I]-x[I])*elementBoundaryLocalOuterNormalsArray[eN*nElementBoundaries_element*nSpace+ ebN*nSpace + I];
		}
	      if (std::abs(dxf) <= boundaryTolerance)
		{
		  for (int I=0; I < nSpace; I++)
		    x[I] += dxf*elementBoundaryLocalOuterNormalsArray[eN*nElementBoundaries_element*nSpace+ ebN*nSpace + I];
		}
	    }//ebN close to boundary
	  bool needToAdd = true;
	  //check points in the element to see if the point is too close
	  if (elementsToTrackedPoints.find(eN) != elementsToTrackedPoints.end())
	    {
	      std::list<double>::iterator eN_it = elementsToTrackedPoints[eN].begin();
	      while (eN_it != elementsToTrackedPoints[eN].end())
		{
		  double dxp = 0.0;
		  for (int I = 0; I < 3; I++)
		    {
		      assert(eN_it != elementsToTrackedPoints[eN].end());
		      const double yI = *eN_it;
		      dxp += (x[I]-yI)*(x[I]-yI);
		      eN_it++;
		    }
		  dxp = std::sqrt(dxp);
		  if (dxp <= neighborTolerance)
		    {
		      needToAdd = false;
		      //mwf debug
		      //std::cout<<"assignTrackedPoints rejecting tracked point ";
		      //for (int I = 0; I < 3; I++)
			//std::cout<<" "<<x[I]<<" ";
		      //std::cout<<std::endl;
		      break;
		    }
		}
	    }
	  else //check to see if the point is too close to one of the vertices
	    {
	      for (int nN_local = 0; nN_local < nNodes_element; nN_local++)
		{
		  const int nN_global = elementNodesArray[eN*nNodes_element+nN_local];
		  double dxp=0.0;
		  for (int I=0; I < 3; I++)
		    {
		      const double yI = nodeArray[nN_global*3+I];
		      dxp += (x[I]-yI)*(x[I]-yI);
		    }
		  dxp = std::sqrt(dxp);
		  if (dxp <= neighborTolerance)
		    {
		      needToAdd = false;
		      break;
		    }
		}
	    }
	  //still need to add
	  if (needToAdd && elementsToTrackedPoints.find(eN) != elementsToTrackedPoints.end())
	    {
	      for (int I=0; I < 3; I++)
		elementsToTrackedPoints[eN].push_back(x[I]);
	    }
	  else if (needToAdd)
	    {
	      std::list<double> tmp; 
	      //add point
	      for (int I=0; I < 3; I++)
		tmp.push_back(x[I]);
	      elementsToTrackedPoints[eN] = tmp;
	      //and the vertices of the element
	      for (int nN=0; nN < nNodes_element; nN++)
		{
		  const int nN_global = elementNodesArray[eN*nNodes_element+nN];
		  for (int I=0; I < 3; I++)
		    elementsToTrackedPoints[eN].push_back(nodeArray[nN_global*3+I]);
		}
	    }
	}//end adding a valid point
    }//k loop over points
}

void createElementSSIPs_1d(const std::list<double>& elementTrackedPoints,
			   std::vector<double>& x_element, 
			   std::vector<double>& w_element)
{
  /***********************************************************************
    assume input is a list of coordinates for a 1d element with SSIPs 
      (and vertices) 
    determine triangulation of element with SSIPs as vertices
    then generate 2nd order gaussian quadrature points and weights on each 
      sub-interval/element of triangulation 
      
   ***********************************************************************/

  unsigned int nPointsIn = elementTrackedPoints.size()/3;
  assert(nPointsIn > 0);

  x_element.clear(); w_element.clear();

  //sort points 
  //if size is a problem could use a fancy comp operator to compare x entries
  std::vector<std::vector<double> > tmpEdgeMesh(nPointsIn);
  std::list<double>::const_iterator it = elementTrackedPoints.begin();
  int ip=0;
  while (it != elementTrackedPoints.end())
    {
      tmpEdgeMesh[ip].resize(3);
      for (int I=0; I < 3; I++)
	{
	  tmpEdgeMesh[ip][I] = *it;
	  it++;
	}
      ip++;
    }
  //I believe default should work
  std::sort(tmpEdgeMesh.begin(),tmpEdgeMesh.end());

  const int nElementsBase = nPointsIn-1;
  const int nSubElementPoints = 2;
  //hardwire 2nd order Gaussian quadrature for now
  double subElementPointsRef[nSubElementPoints] = {(sqrt(3.0)-1.0)/(2.0*sqrt(3.0)),
						   (sqrt(3.0)+1.0)/(2.0*sqrt(3.0))};
  double subElementWeightsRef[nSubElementPoints]= {0.5,0.5};

  const int nQuadraturePointsNew = nElementsBase*nSubElementPoints;

  double d[3] = {0.0,0.0,0.0}; 

  for (int eN_local=0; eN_local < nElementsBase; eN_local++)
    {
      double volume = 0.0;
      for (int I=0; I < 3; I++)
	{
	  d[I] = tmpEdgeMesh[eN_local+1][I]-tmpEdgeMesh[eN_local][I];
	  //mwf debug
	  //if (I==0)
	  // assert(d[I] > 0.0);
	  volume += d[I]*d[I];
	}
      volume = std::sqrt(volume);
      for (int k=0; k < nSubElementPoints; k++)
	{
	  for (int I=0; I < 3; I++)
	    {
	      const double xk = tmpEdgeMesh[eN_local][I]*(1.0-subElementPointsRef[k]) + tmpEdgeMesh[eN_local+1][I]*subElementPointsRef[k];
	      x_element.push_back(xk);
	    }
	  w_element.push_back(volume*subElementWeightsRef[k]);
	}//sub-quadrature points
    }//eN_local loop over element triangulation
  assert(x_element.size() == unsigned(nQuadraturePointsNew*3));
  assert(w_element.size() == unsigned(nQuadraturePointsNew));
  for (int i=0; i < nQuadraturePointsNew-1; i++)
    {
      assert(x_element[(i+1)*3] > x_element[i*3]);
    }
}

//----------------------------------------------------------------------
// Postprocess just particle trajectory results to get concentration
// information using convolution-based particle tracking (CBPT) approach
//
extern "C"
double integratePiecewiseLinearMassSource(int nknots,
					 const double * t_vals,
					 const double * m_vals,
					 double t_in,
					 double t_out,
					 double tau_out,
					 double tau_in,
					 double decay) //must be nonpositive
{
  //TODO replace with faster search
  double massint = 0.0; //return value
  int i = 0;
  assert (tau_in > tau_out);
  while (i < nknots-1 && t_vals[i] < tau_out) i++;
  if (i > 0) i--; //i is start of interval containing tau_out
  assert(i == nknots - 2 || t_vals[i] <= tau_out);

  int j = i+1;
  while (j < nknots-1 && t_vals[j] < tau_in) j++;
  //j should be end of interval containing tau_in
  assert(j == nknots-1 || t_vals[j] >= tau_in); 

  int nsegments = j-i;
  double decay_inv_abs = 0.0;
  if (decay < 0.0)
    decay_inv_abs = 1.0/fabs(decay);

  double tau_m = tau_out;  
  for (int k = 0; k  < nsegments; k++)
    {
      //slope over segment
      double beta = (m_vals[i+k+1] - m_vals[i+k])/(t_vals[i+k+1]-t_vals[i+k]);
      double mass_m= m_vals[i+k] + beta*(tau_m - t_vals[i+k]);

      double tau_p = fmin(tau_in,t_vals[i+k+1]);
      double mass_p= m_vals[i+k] + beta*(tau_p - t_vals[i+k]);
      //corresponding values of t for evaluating exponentials
      double t_m = (tau_p - tau_out) + t_in; // --> t_out
      double t_p = (tau_m - tau_out) + t_in; // --> t_in
      if (decay < 0.0) //equation 16
	{
	  massint += (mass_p - beta*decay_inv_abs)*exp(decay*t_p)*decay_inv_abs - 
	    (mass_m - beta*decay_inv_abs)*exp(decay*t_m)*decay_inv_abs;
	}
      else
	{
	  massint += (tau_p-tau_m)*(mass_m+mass_p)*0.5;
	}
      tau_m = tau_p;
    }
  return massint;
}
extern "C"
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
				  double *c_element) //element concentrations
{
  const double particle_weight = 1.0/float(nParticles_global);
  //add as argument
  const double dt_tol = 1.0e-8;
  for (int k = 0; k < nParticles_global; k++)
    {
      int flag_k = particleFlags[k];
      int i = traj_offsets[k]; 
      double t_in = t_traj[i]; 
      int eN = elem_traj[i]; //current element
      int i_in = i;
      double decay = 0.0;
      while (i < traj_offsets[k+1]) //walk through trajectory for particle k
	{
	  //physical coefficients for this element
	  decay += decay_coef_element[eN*nParticleFlags + flag_k];
	  double retardation = retardation_factor_element[eN*nParticleFlags + flag_k];
	  assert (fabs(retardation) > 0.0);
	  while (i < traj_offsets[k+1] && elem_traj[i] == eN) //walk until through element
            i++;
          int i_out = i-1; if (i_out < i_in) i_out = i_in; //mwf check this
	  //adjust travel times due to retardation?
	  double dt_cons  = t_traj[i_out]-t_traj[i_in];//travel time for conservative simulation
	  double t_out    = t_traj[i_out];
	  t_out = t_in + dt_cons*retardation; 

	  if (fabs(t_out-t_in) > dt_tol)
	    {
	      
	      double tau_out = fmax(0.,tau-t_out), tau_in = fmax(0.0,tau-t_in);
	      double convolution_term = 0.0;
	      if (tau_out < tau_in) //evaluate mass source
		convolution_term = integratePiecewiseLinearMassSource(nMassSourceKnots,massSource_t,massSource_m,t_in,t_out,tau_out,tau_in,decay);
	      c_element[eN] += convolution_term*particle_weight;
	    }
	  //get ready for next element
	  assert(eN == elem_traj[i_out] || i == traj_offsets[k+1]);
          assert(eN != elem_traj[i]);
          
          eN = elem_traj[i];
	  t_in  = t_out;
          i_in  = i_out;
          
	}// end walk through trajectory for i
    }//end particle loop
}
extern "C"
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
					       double *c_element) //element concentrations
{
  const double particle_weight = 1.0/float(nParticles_global);
  //add as argument
  const double dt_tol = 1.0e-8;
  for (int k = 0; k < nParticles_global; k++)
    {
      int flag_k = particle_flags[k];
      int i = traj_offsets[k]; 
      double t_in = t_traj[i]; 
      int eN = elem_traj[i]; //current element
      int i_in = i;
      double decay = 0.0;
      while (i < traj_offsets[k+1]) //walk through trajectory for particle k
	{
	  //physical coefficients for this element
	  const int material_id = material_types_element[eN];
	  decay += decay_coef_types[material_id*nParticleFlags + flag_k];
	  double retardation = retardation_factor_types[material_id*nParticleFlags + flag_k];
	  assert (fabs(retardation) > 0.0);
	  while (i < traj_offsets[k+1] && elem_traj[i] == eN) //walk until through element
            i++;
          int i_out = i-1; if (i_out < i_in) i_out = i_in; //mwf check this
	  //adjust travel times due to retardation?
	  double dt_cons  = t_traj[i_out]-t_traj[i_in];//travel time for conservative simulation
	  double t_out    = t_traj[i_out];
	  t_out = t_in + dt_cons*retardation; 
	  //mwf debug
	  //std::cout<<"accumSourceMat i= "<<i<<" eN= "<<eN<<" matid= "<<material_id<<" flag_k= "<<flag_k<<" nParticleFlags= "<<nParticleFlags<<" decay = "<<decay<<" R= "
	  //	   <<" t_in= "<<t_in<<" t_out= "<<t_out<<retardation<<std::endl;
	  if (fabs(t_out-t_in) > dt_tol)
	    {
	      
	      double tau_out = fmax(0.,tau-t_out), tau_in = fmax(0.0,tau-t_in);
	      double convolution_term = 0.0;
	      if (tau_out < tau_in) //evaluate mass source
		convolution_term = integratePiecewiseLinearMassSource(nMassSourceKnots,massSource_t,massSource_m,t_in,t_out,tau_out,tau_in,decay);
	      //mwf debug
	      //std::cout<<"\t"<<" tau_out= "<<tau_out<<" tau_in= "<<tau_in<<" conv_term= "<<convolution_term<<std::endl;
	      c_element[eN] += convolution_term*particle_weight;
	    }
	  //get ready for next element
	  assert(eN == elem_traj[i_out] || i == traj_offsets[k+1]);
          assert(eN != elem_traj[i]);
          
          eN = elem_traj[i];
	  t_in  = t_out;
          i_in  = i_out;
	  assert(eN >= 0 && eN < nElements_global);
	}// end walk through trajectory for i
    }//end particle loop
}


extern "C"
double volume123(int nSpace, //space dimension
		 int nNodes_element, //number of nodes in this element
		 const int* elementNodes_element, //element--> global node table for this element 
		 const double* nodeArray) //coordinates (3d)
{
  //use Pearce's approach to compute volume for simplex, quad, prism, or hex
  //base 1
  const int KVB2[2][2][3] = {{{1,2,3},{1,3,4}},  
			     {{1,2,3},{0,0,0}}};
  //base 1
  const int KVB3[3][6][4] = {{{5,1,2,4}, {6,5,2,4}, {8,5,6,4}, {7,2,3,4}, {8,6,7,4}, {6,7,4,2}},
			     {{4,1,2,3}, {5,4,2,3}, {6,4,5,3}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}},
			     {{4,1,2,3}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}}};
  //base 1
  const int KKOFF[4][3] = {{2,3,4},{1,3,4},{1,2,4},{1,2,3}}; 
  const int ibase = 1;
  double XX[4] = {0.,0.,0.,0.}, YY[4] = {0.,0.,0.,0.}, ZZ[4] = {0.,0.,0.,0.};
  double volume=0.0;
  if (nSpace == 1)
    {
      assert(nNodes_element == 2);
      volume = fabs(nodeArray[elementNodes_element[1]*3+0]-
		    nodeArray[elementNodes_element[0]*3+0]);
      
      return volume;
    }//1d
  if (nSpace == 2)
    {
      //default is triangle
      int ID = 1, NA=1; //which slice of KVB2 and how many sub triangles
      if (nNodes_element == 4)
	{
	  ID = 0; NA=2; //quad
	}
      for (int IA =0 ; IA < NA; IA++) //loop through sub-triangles
	{
	  for (int I=0; I < 3; I++) //build triangle
	    {
	      const int II = KVB2[ID][IA][I]-ibase;
	      assert(0 <= II && II < nNodes_element);
	      const int nN = elementNodes_element[II]; 
	      XX[I] = nodeArray[nN*3 + 0];
	      YY[I] = nodeArray[nN*3 + 1];
	    }
	  //determinant
	  const double X12=XX[0]-XX[1];
	  const double X23=XX[1]-XX[2];
          const double X31=XX[2]-XX[0]; 
          const double Y12=YY[0]-YY[1];
          const double Y23=YY[1]-YY[2];
          const double Y31=YY[2]-YY[0];
	  volume += XX[0]*Y23 + XX[1]*Y31 + XX[2]*Y12;
	}
      volume = fabs(volume)/double(NA);
      return volume;
    }//2d
  if (nSpace == 3)
    {

      int ID = 2, NV=1; //which slice of KVB3 and how many sub tets
      if (nNodes_element == 6) //triangular prism
	{
	  ID=1; NV = 3;
	}
      else if (nNodes_element == 8) //hexahedron
	{
	  ID = 0; NV = 6;
	}
      for (int IV = 0; IV < NV; IV++) //loop through sub-tets
	{
	  for (int I = 0; I < 4; I++) //build tet
	    {
	      const int II = KVB3[ID][IV][I]-ibase;
	      const int nN = elementNodes_element[II];
	      XX[I] = nodeArray[nN*3+0]; YY[I] = nodeArray[nN*3+1]; ZZ[I] = nodeArray[nN*3+2];
	    }
	  for (int KK=0; KK < 4; KK++)
	    {
	      const int K1 = KKOFF[KK][0]-ibase, K2 = KKOFF[KK][1]-ibase, K3 = KKOFF[KK][2]-ibase;
              const double tmp = pow(-1.0,KK)*(XX[K1]*YY[K2]*ZZ[K3]+
                                               YY[K1]*ZZ[K2]*XX[K3]+ZZ[K1]*XX[K2]*YY[K3]-
                                               XX[K3]*YY[K2]*ZZ[K1]-YY[K3]*ZZ[K2]*XX[K1]-
                                               ZZ[K3]*XX[K2]*YY[K1]);
	      volume += tmp;
	    }//KK
	}//IV
      volume = fabs(volume)/float(NV);
      return volume;
    }
}
extern "C"
void nodalProjection123(int calculateNodeStarVolume,
			int nElements_global,
			int nNodes_global,
			const int* elementNodes_offset,
			const int* elementNodesArray,
			const double* elementVolumesArray,
			const double* c_ele,
			double * nodeStarVolume, 
			double * c_nodal)
{
  //mass-lumped L2 projection of element-wise constant concentration
  //loop through all the elements and accumulate that elements contribution to 
  //all of its nodes
  memset(c_nodal,0,sizeof(double)*nNodes_global);
  if (calculateNodeStarVolume > 0)
    {
      memset(nodeStarVolume,0,sizeof(double)*nNodes_global);
      for (int eN=0; eN < nElements_global; eN++)
	{
	  const int nNodes_element = elementNodes_offset[eN+1] - elementNodes_offset[eN];
	  const double weight = elementVolumesArray[eN]/double(nNodes_element);
	  for (int i=elementNodes_offset[eN]; i < elementNodes_offset[eN+1]; i++)
	    {
	      const int nN = elementNodesArray[i];
	      nodeStarVolume[nN] += weight;
	    }
	}
    }
  
  for (int eN=0; eN < nElements_global; eN++)
    {
      const int nNodes_element = elementNodes_offset[eN+1] - elementNodes_offset[eN];
      const double weight = elementVolumesArray[eN]/double(nNodes_element);
      for (int i=elementNodes_offset[eN]; i < elementNodes_offset[eN+1]; i++)
	{
	  const int nN = elementNodesArray[i];
	  c_nodal[nN] += weight*c_ele[eN];
	}
    }
  for (int nN = 0; nN < nNodes_global; nN++)
    {
      c_nodal[nN] /= nodeStarVolume[nN];
    }
}
}//namespace ELLAM
