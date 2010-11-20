#include "ellam.h"
#include "tracking.h"
#include <cmath>
#include <iostream>
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
//			int offset_u, int stride_u, 
//			double* globalResidual)
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
       5. accumulate m^{n+1}_k*W_k*w_i(x^{n+1}_{k}) to global residual I(i) 
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
					 const double* dm,
					 const double* w,
					 const double* v,
					 const double* dV,
					 const double* rhs,
					 double* elementResidual,
					 double* theta,
					 double* slumpedMassMatrix)
{
  int eN,i,j,k,I,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  int stencilWidth=5;
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
      
      /*calculate local mass matrix, approximation isn't exactly right
	for now, since assumes matrix same on both elements*/
      for (i=0; i < nDOF_test_element; i++)
	for (j=0; j < nDOF_trial_element; j++)
	  {
	    slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = 0.0;
	    for (k=0; k < nQuadraturePoints_element; k++)
	      {
		slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] +=  
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
      a = slumpedMassMatrix[eN*nDOF_test_X_trial_element + 0*nDOF_trial_element + 0];
      b = slumpedMassMatrix[eN*nDOF_test_X_trial_element + 0*nDOF_trial_element + 1];

      tmp1 = dr0/(drm1+dr0+drp1+1.0e-12);
      tmp2 = fmin(drm1,dr0)/(drm1+dr0+1.0e-12);
      tmp3 = fmin(drp1,dr0)/(drp1+dr0+1.0e-12);
      tmp4 = fmin(fmin(tmp1,tmp2),tmp3);
      theta[eN] = b - tmp4*(a+b);
      theta[eN] = fmax(theta[eN],0.0);

      /*update mass matrix and element residual*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] += theta[eN];
	  elementResidual[eN*nDOF_test_element + i] += theta[eN]*u_dof[l2g[eN*nDOF_trial_element+i]];

	  for (j=0; j < i; j++)
	    {
	      slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] -= theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	  for (j=i+1; j < nDOF_trial_element; j++)
	    {
	      slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] -= theta[eN];
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
					 const int* l2g,
					 const int* elementNeighborsArray,
					 const double* u_dof,
					 const double* dm,
					 const double* w,
					 const double* v,
					 const double* dV,
					 const double* rhs,
					 double* elementResidual,
					 double* theta,
					 double* slumpedMassMatrix)
{
  /*
    try stupid approach where just slump local element mass matrix
   */
  int eN,i,j,k,I,nDOF_test_X_trial_element=nDOF_trial_element*nDOF_test_element;
  
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
      /*calculate local mass matrix */
      for (i=0; i < nDOF_test_element; i++)
	for (j=0; j < nDOF_trial_element; j++)
	  {
	    slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] = 0.0;
	    for (k=0; k < nQuadraturePoints_element; k++)
	      {
		slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] +=  
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
	  a = fmax(a,slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i]);
	  for (j=0; j < nDOF_trial_element; j++)
	    {
	      b = fmin(b,slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j]);
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
      /*update mass matrix and element residual*/
      for (i=0; i < nDOF_test_element; i++)
	{
	  slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + i] += (nDOF_trial_element-1)*theta[eN];
	  elementResidual[eN*nDOF_test_element + i] += (nDOF_trial_element-1)*theta[eN]*u_dof[l2g[eN*nDOF_trial_element+i]];

	  for (j=0; j < i; j++)
	    {
	      slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] -= theta[eN];
	      elementResidual[eN*nDOF_test_element + i] -= theta[eN]*u_dof[l2g[eN*nDOF_trial_element+j]];
	    }
	  for (j=i+1; j < nDOF_trial_element; j++)
	    {
	      slumpedMassMatrix[eN*nDOF_test_X_trial_element + i*nDOF_trial_element + j] -= theta[eN];
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
