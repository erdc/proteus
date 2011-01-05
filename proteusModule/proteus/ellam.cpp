#include "ellam.h"
#include "tracking.h"
#include <cmath>
#include <iostream>
#include <cassert>
#include <list>
#include <map>

#include PROTEUS_LAPACK_H
/***********************************************************************

  TODO:

   switch dimension of dV and x_track to just be the number of global
   integration points so that there can be different numbers per element

   decide if ignoring outflow mass is a better idea
**********************************************************************/

extern "C" 
void updateOldMass_weak(int nSpace,            
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
	  if (flag_track[eN_k] >= -1)
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
		  //globalResidual[offset_u+stride_u*u_l2g[eN_track_i]] -= m_old*w[i]*weight; 
		  //mwf debug
		  //std::cout<<"LADRellam oldMass globalResid["<<offset_u+stride_u*u_l2g[eN_track_i]<<"]= "<< globalResidual[offset_u+stride_u*u_l2g[eN_track_i]] 
		  //   <<" oldmass = "<<m_old*w[i]*weight <<std::endl;
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element
    }//nElements loop
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
      if (flag_track[k] >= -1)
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
void updateExteriorOutflowBoundaryFlux(double dtnp1,          //full time step size
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
	      //mwf debug
	      //std::cout<<"LADR2Dellam boundaryResidual globalResid["<<offset_u+stride_u*u_l2g[eN_i]<<"]= "<< globalResidual[offset_u+stride_u*u_l2g[eN_i]] 
	      //   <<" elementResidual = "<<elementResidual_u[i] <<std::endl;
	    }//i
	}//kb
    }//ebNE
}
extern "C"
void updateExteriorOutflowBoundaryFluxInGlobalResidual(double dtnp1,          //full time step size
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
	      //mwf debug
	      //std::cout<<"LADR2Dellam boundaryResidual globalResid["<<offset_u+stride_u*u_l2g[eN_i]<<"]= "<< globalResidual[offset_u+stride_u*u_l2g[eN_i]] 
	      //   <<" elementResidual = "<<elementResidual_u[i] <<std::endl;
	    }//i
	}//kb
    }//ebNE
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
void accumulateInflowFluxInGlobalResidual(int nSpace,
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
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element boundary
    }//element boundaries
  
}


/**********************************************************************
  go through inflow boundary points that have been tracked forward
  and accumulate the inflow flux into the element residual for the test functions
  with non-zero support in those areas
 **********************************************************************/
extern "C" 
void accumulateInflowFlux(int nSpace,
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
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element boundary
    }//element boundaries
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
      
      double max_delta_IJ = 0.0; double min_delta_IJ = 0.0;
      assert(nodeStarOffsets[I0+1] - nodeStarOffsets[I0] > 1);
      int offset_max_delta_IJ  = nodeStarOffsets[I0]; int offset_min_delta_IJ = nodeStarOffsets[I0]+1; 
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
	  if (delta_IJ > max_delta_IJ)
	    {
	      offset_max_delta_IJ = offset;
	      max_delta_IJ = delta_IJ;
	    }
	  if (delta_IJ < min_delta_IJ)
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
	  if (delta_Ip1J > max_delta_Ip1J)
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

void generateArraysForSSIPs(int nElements_global,
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
  ELLAM::assignTrackedPointsToElements(nElements_global,
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
	  unsigned int nPoints_eN = 0;
	  if (elementsToTrackedPoints.find(eN) != elementsToTrackedPoints.end())
	    {
	      ELLAM::createElementSSIPs_1d(elementsToTrackedPoints[eN],
					   x_element,
					   w_element);
	      nPoints_eN = w_element.size();
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
      std::cout<<"generateArraysForSSIPs nSpace= "<<nSpace<<" not implemented! "<<std::endl;
      assert(0);
    }

}
namespace ELLAM
{
void assignTrackedPointsToElements(int nElements_global,
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
}//namespace ELLAM
