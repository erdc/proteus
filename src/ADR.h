#ifndef ADR_H
#define ADR_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

namespace proteus
{
  class cppADR_base
  {
  public:
    virtual ~cppADR_base(){}
    virtual void calculateResidual(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   int* mesh_l2g,
				   double* dV_ref,
				   double* u_trial_ref,
				   double* u_grad_trial_ref,
				   double* u_test_ref,
				   double* u_grad_test_ref,
				   //element boundary
				   double* mesh_trial_trace_ref,
				   double* mesh_grad_trial_trace_ref,
				   double* dS_ref,
				   double* u_trial_trace_ref,
				   double* u_grad_trial_trace_ref,
				   double* u_test_trace_ref,
				   double* u_grad_test_trace_ref,
				   double* normal_ref,
				   double* boundaryJac_ref,
				   //physics
				   int nElements_global,
				   int* u_l2g, 
				   double* u_dof,
				   int* sd_rptr,
				   int* sd_colind,
				   double* q_a,
				   double* q_r,
				   int offset_u, int stride_u, 
				   double* globalResidual,			   
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   double* ebqe_a,
				   int* isDOFBoundary_u,
				   double* ebqe_bc_u_ext,
				   int* isDiffusiveFluxBoundary_u,
				   double* ebqe_bc_flux_u_ext,
				   double* ebqe_penalty_ext)=0;
    virtual void calculateJacobian(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   int* mesh_l2g,
				   double* dV_ref,
				   double* u_trial_ref,
				   double* u_grad_trial_ref,
				   double* u_test_ref,
				   double* u_grad_test_ref,
				   //element boundary
				   double* mesh_trial_trace_ref,
				   double* mesh_grad_trial_trace_ref,
				   double* dS_ref,
				   double* u_trial_trace_ref,
				   double* u_grad_trial_trace_ref,
				   double* u_test_trace_ref,
				   double* u_grad_test_trace_ref,
				   double* normal_ref,
				   double* boundaryJac_ref,
				   //physics
				   int nElements_global,
				   int* u_l2g,
				   double* u_dof,
				   int* sd_rowptr,
				   int* sd_colind,
				   double* q_a,
				   double* q_r,
				   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				   double* globalJacobian,
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   double* ebqe_a,
				   int* isDOFBoundary_u,
				   double* ebqe_bc_u_ext,
				   int* isDiffusiveFluxBoundary_u,
				   double* ebqe_bc_flux_u_ext,
				   int* csrColumnOffsets_eb_u_u,
				   double* ebqe_penalty_ext)=0;
  };
  
  template<class CompKernelType,
    int nSpace,
    int nQuadraturePoints_element,
    int nDOF_mesh_trial_element,
    int nDOF_trial_element,
    int nDOF_test_element,
    int nQuadraturePoints_elementBoundary>
    class cppADR : public cppADR_base
    {
    public:
      const int nDOF_test_X_trial_element;
      CompKernelType ck;
    cppADR():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
	ck()
	  {}
      /* inline
	 void evaluateCoefficients()
	 {
	 
	 }
      */
      
      inline
	void exteriorNumericalDiffusiveFlux(int* rowptr,
					  int* colind,
					  const int& isDOFBoundary,
					  const int& isDiffusiveFluxBoundary,
					  const double n[nSpace],
					  double* bc_a,
					  const double& bc_u,
					  const double& bc_flux,
					  double* a,
					  const double grad_potential[nSpace],
					  const double& u,
					  const double& penalty,
					  double& flux)
    {
      double diffusiveVelocityComponent_I,penaltyFlux,max_a;
      if(isDiffusiveFluxBoundary == 1)
	{
	  flux = bc_flux;
	}
      else if(isDOFBoundary == 1)
	{
	  flux = 0.0;
	  max_a=0.0;
	  for(int I=0;I<nSpace;I++)
	    {
	      diffusiveVelocityComponent_I=0.0;
	      for(int m=rowptr[I];m<rowptr[I+1];m++)
		{
		  diffusiveVelocityComponent_I -= a[m]*grad_potential[colind[m]];
		  max_a = fmax(max_a,a[m]);
		}
	      flux+= diffusiveVelocityComponent_I*n[I];
	    }
	  penaltyFlux = max_a*penalty*(u-bc_u);
	  std::cout<<"max_a "<<max_a<<" penalty "<<penalty<<" u "<<u<<" bc_u "<<bc_u<<std::endl;
	  std::cout<<"flux "<<flux<<" penalty flux "<<penaltyFlux<<std::endl;
	  flux += penaltyFlux;
	}
      else
	{
	  std::cerr<<"warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl;
	  flux = 0.0;
	}
    }


    inline
    double ExteriorNumericalDiffusiveFluxJacobian(int* rowptr,
						  int* colind,
						  const int& isDOFBoundary,
						  const int& isDiffusiveFluxBoundary,
						  const double n[nSpace],
						  double* a,
						  const double& v,
						  const double grad_v[nSpace],
						  const double& penalty)
    {
      double dvel_I,tmp=0.0,max_a=0.0;
      if(isDiffusiveFluxBoundary==0 && isDOFBoundary==1)
	{
	  for(int I=0;I<nSpace;I++)
	    {
	      dvel_I=0.0;
	      for(int m=rowptr[I];m<rowptr[I+1];m++)
		{
		  dvel_I -= a[m]*grad_v[colind[m]];
		  max_a = fmax(max_a,a[m]);
		}
	      tmp += dvel_I*n[I];
	    }
	  tmp +=max_a*penalty*v;
	}
      return tmp;
    }

    inline void calculateElementResidual(//element
					 double* mesh_trial_ref,
					 double* mesh_grad_trial_ref,
					 double* mesh_dof,
					 int* mesh_l2g,
					 double* dV_ref,
					 double* u_trial_ref,
					 double* u_grad_trial_ref,
					 double* u_test_ref,
					 double* u_grad_test_ref,
					 //element boundary
					 double* mesh_trial_trace_ref,
					 double* mesh_grad_trial_trace_ref,
					 double* dS_ref,
					 double* u_trial_trace_ref,
					 double* u_grad_trial_trace_ref,
					 double* u_test_trace_ref,
					 double* u_grad_test_trace_ref,
					 double* normal_ref,
					 double* boundaryJac_ref,
					 //physics
					 int nElements_global,
					 int* u_l2g, 
					 double* u_dof,
					 int* sd_rowptr,
					 int* sd_colind,
					 double* q_a,
					 double* q_r,
					 int offset_u, int stride_u, 
					 double* elementResidual_u,			   
					 int nExteriorElementBoundaries_global,
					 int* exteriorElementBoundariesArray,
					 int* elementBoundaryElementsArray,
					 int* elementBoundaryLocalElementBoundariesArray,
					 double* element_u,
					 int eN)
    {
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_u[i]=0.0;
	}//i
      //loop over quadrature points and compute integrands
      for  (int k=0;k<nQuadraturePoints_element;k++)
	{
	  //compute indeces and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k;
	  register double u=0.0,grad_u[nSpace],
	    *a=NULL,
	    r=0.0,
	    jac[nSpace*nSpace],
	    jacDet,
	    jacInv[nSpace*nSpace],
	    u_grad_trial[nDOF_trial_element*nSpace],
	    u_test_dV[nDOF_trial_element],
	    u_grad_test_dV[nDOF_test_element*nSpace],
	    dV,x,y,z,
	    G[nSpace*nSpace],G_dd_G,tr_G;
	  //
	  //compute solution and gradients at quadrature points
	  //
	  ck.calculateMapping_element(eN,
				      k,
				      mesh_dof,
				      mesh_l2g,
				      mesh_trial_ref,
				      mesh_grad_trial_ref,
				      jac,
				      jacDet,
				      jacInv,
				      x,y,z);
	  //get the physical integration weight
	  dV = fabs(jacDet)*dV_ref[k];
	  //get the metric tensor and friends
	  ck.calculateG(jacInv,G,G_dd_G,tr_G);
	  //get the trial function gradients
	  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	  //get the solution
	  ck.valFromElementDOF(element_u,&u_trial_ref[k*nDOF_trial_element],u);
	  //get the solution gradients
	  ck.gradFromElementDOF(element_u,u_grad_trial,grad_u);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
	      for (int I=0;I<nSpace;I++)
		{
		  u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	    }
	  
	  //
	  //calculate pde coefficients at quadrature points
	  //
	  //evaluateCoefficients();
	  //just set from pre-evaluated quadrature point values for now
	  a = &q_a[eN_k*sd_rowptr[nSpace]];
	  r = q_r[eN_k];
	  // 
	  //update element residual 
	  // 
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int  i_nSpace=i*nSpace;
	      
	      elementResidual_u[i] += ck.Diffusion_weak(sd_rowptr,sd_colind,a,grad_u,&u_grad_test_dV[i_nSpace]) + 
		ck.Reaction_weak(r,u_test_dV[i]);
	    }//i
	}
    }

    void calculateResidual(//element
			   double* mesh_trial_ref,
			   double* mesh_grad_trial_ref,
			   double* mesh_dof,
			   int* mesh_l2g,
			   double* dV_ref,
			   double* u_trial_ref,
			   double* u_grad_trial_ref,
			   double* u_test_ref,
			   double* u_grad_test_ref,
			   //element boundary
			   double* mesh_trial_trace_ref,
			   double* mesh_grad_trial_trace_ref,
			   double* dS_ref,
			   double* u_trial_trace_ref,
			   double* u_grad_trial_trace_ref,
			   double* u_test_trace_ref,
			   double* u_grad_test_trace_ref,
			   double* normal_ref,
			   double* boundaryJac_ref,
			   //physics
			   int nElements_global,
			   int* u_l2g, 
			   double* u_dof,
			   int* sd_rowptr,
			   int* sd_colind,
			   double* q_a,
			   double* q_r,
			   int offset_u, int stride_u, 
			   double* globalResidual,			   
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
			   double* ebqe_a,
			   int* isDOFBoundary_u,
			   double* ebqe_bc_u_ext,
			   int* isDiffusiveFluxBoundary_u,
			   double* ebqe_bc_flux_u_ext,
			   double* ebqe_penalty_ext)
    {
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      //eN is the element index
      //eN_k is the quadrature point index for a scalar
      //eN_k_nSpace is the quadrature point index for a vector
      //eN_i is the element test function index
      //eN_j is the element trial function index
      //eN_k_j is the quadrature point index for a trial function
      //eN_k_i is the quadrature point index for a trial function
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for element residual and initialize
	  register double elementResidual_u[nDOF_test_element],element_u[nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i=eN*nDOF_test_element+i;
	      element_u[i] = u_dof[u_l2g[eN_i]];
	    }//i
	  calculateElementResidual(mesh_trial_ref,
				   mesh_grad_trial_ref,
				   mesh_dof,
				   mesh_l2g,
				   dV_ref,
				   u_trial_ref,
				   u_grad_trial_ref,
				   u_test_ref,
				   u_grad_test_ref,
				   mesh_trial_trace_ref,
				   mesh_grad_trial_trace_ref,
				   dS_ref,
				   u_trial_trace_ref,
				   u_grad_trial_trace_ref,
				   u_test_trace_ref,
				   u_grad_test_trace_ref,
				   normal_ref,
				   boundaryJac_ref,
				   nElements_global,
				   u_l2g, 
				   u_dof,
				   sd_rowptr,
				   sd_colind,
				   q_a,
				   q_r,
				   offset_u,stride_u, 
				   elementResidual_u,			   
				   nExteriorElementBoundaries_global,
				   exteriorElementBoundariesArray,
				   elementBoundaryElementsArray,
				   elementBoundaryLocalElementBoundariesArray,
				   element_u,
				   eN);
	  //
	  //load element into global residual and save element residual
	  //
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_i=eN*nDOF_test_element+i;
          
	      globalResidual[offset_u+stride_u*u_l2g[eN_i]]+=elementResidual_u[i];
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE], 
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_u[i]=0.0;
	    }
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	      register double u_ext=0.0,
		grad_u_ext[nSpace],
		m_ext=0.0,
		dm_ext=0.0,
		*a_ext,
		/* *da_exxt, */
		/* f_ext[nSpace], */
		/* df_ext[nSpace], */
		r_ext=0.0,
		/* dr_ext=0.0, */
		flux_diff_ext=0.0,
		bc_u_ext=0.0,
		//bc_grad_u_ext[nSpace],
		/* bc_m_ext=0.0, */
		/* bc_dm_ext=0.0, */
		/* bc_f_ext[nSpace], */
		/* bc_df_ext[nSpace], */
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		u_grad_trial_trace[nDOF_trial_element*nSpace],
		u_grad_test_dS[nDOF_trial_element*nSpace],
		normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		//
		G[nSpace*nSpace],G_dd_G,tr_G;
	      // 
	      //calculate the solution and gradients at quadrature points 
	      // 
	      //compute information about mapping from reference element to physical element
	      ck.calculateMapping_elementBoundary(eN,
						  ebN_local,
						  kb,
						  ebN_local_kb,
						  mesh_dof,
						  mesh_l2g,
						  mesh_trial_trace_ref,
						  mesh_grad_trial_trace_ref,
						  boundaryJac_ref,
						  jac_ext,
						  jacDet_ext,
						  jacInv_ext,
						  boundaryJac,
						  metricTensor,
						  metricTensorDetSqrt,
						  normal_ref,
						  normal,
						  x_ext,y_ext,z_ext);
	      /* ck.calculateMappingVelocity_elementBoundary(eN, */
	      /* 						  ebN_local, */
	      /* 						  kb, */
	      /* 						  ebN_local_kb, */
	      /* 						  mesh_velocity_dof, */
	      /* 						  mesh_l2g, */
	      /* 						  mesh_trial_trace_ref, */
	      /* 						  xt_ext,yt_ext,zt_ext, */
	      /* 						  normal, */
	      /* 						  boundaryJac, */
	      /* 						  metricTensor, */
	      /* 						  integralScaling); */
	      /*dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];*/
	      dS = metricTensorDetSqrt*dS_ref[kb];
	      //get the metric tensor
	      //cek todo use symmetry
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  for (int I=0;I<nSpace;I++)
		    u_grad_test_dS[j*nSpace+I] = u_grad_trial_trace[j*nSpace+I]*dS;//cek hack, using trial
		}
	      //
	      //load the boundary values
	      //
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	      //
	      // 
	      //calculate the pde coefficients using the solution and the boundary values for the solution 
	      // 
	      a_ext = &ebqe_a[ebNE_kb*sd_rowptr[nSpace]];
	      /* evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace], */
	      /* 			   u_ext, */
	      /* 			   //VRANS */
	      /* 			   porosity_ext, */
	      /* 			   // */
	      /* 			   m_ext, */
	      /* 			   dm_ext, */
	      /* 			   f_ext, */
	      /* 			   df_ext); */
	      /* evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace], */
	      /* 			   bc_u_ext, */
	      /* 			   //VRANS */
	      /* 			   porosity_ext, */
	      /* 			   // */
	      /* 			   bc_m_ext, */
	      /* 			   bc_dm_ext, */
	      /* 			   bc_f_ext, */
	      /* 			   bc_df_ext);     */
	      //
	      //moving mesh
	      //
	      /* double velocity_ext[nSpace]; */
	      /* double mesh_velocity[3]; */
	      /* mesh_velocity[0] = xt_ext; */
	      /* mesh_velocity[1] = yt_ext; */
	      /* mesh_velocity[2] = zt_ext; */
	      /* for (int I=0;I<nSpace;I++) */
	      /* 	velocity_ext[I] = ebqe_velocity_ext[ebNE_kb_nSpace+0] - MOVING_DOMAIN*mesh_velocity[I]; */
	      // 
	      //calculate the numerical fluxes 
	      // 
	      exteriorNumericalDiffusiveFlux(sd_rowptr,
					     sd_colind,
					     isDOFBoundary_u[ebNE_kb],
					     isDiffusiveFluxBoundary_u[ebNE_kb],
					     normal,
					     a_ext,
					     bc_u_ext,
					     ebqe_bc_flux_u_ext[ebNE_kb],
					     a_ext,
					     grad_u_ext,
					     u_ext,
					     ebqe_penalty_ext[ebNE_kb],
					     flux_diff_ext);
	      std::cout<<"ebNE "<<ebNE<<"kb "<<kb<<"flux_diff_ext"<<flux_diff_ext<<std::endl;
	      /* exteriorNumericalAdvectiveFlux(isDOFBoundary_u[ebNE_kb], */
	      /* 				     isFluxBoundary_u[ebNE_kb], */
	      /* 				     normal, */
	      /* 				     bc_u_ext, */
	      /* 				     ebqe_bc_flux_u_ext[ebNE_kb], */
	      /* 				     u_ext, */
	      /* 				     df_ext, */
	      /* 				     flux_ext); */
	      //ebqe_flux[ebNE_kb] = flux_ext;
	      //
	      //update residuals
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  //int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
		  double eb_adjoint_sigma=0.0;
		  elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_diff_ext,u_test_dS[i]);
 /* + */
 /* 		    ck.ExteriorElementBoundaryDiffusionAdjoint(isDOFBoundary_u[ebNE_kb], */
 /* 							       isDiffusiveFluxBoundary_u[ebNE_kb], */
 /* 							       eb_adjoint_sigma, */
 /* 							       u_ext, */
 /* 							       bc_u_ext, */
 /* 							       normal, */
 /* 							       sd_rowptr, */
 /* 							       sd_colind, */
 /* 							       a_ext, */
 /* 							       &u_grad_test_dS[i*nSpace]); */
		}//i
	    }//kb
	  //
	  //update the element and global residual storage
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
	    }//i
	}//ebNE
    }

    inline void calculateElementJacobian(//element
					 double* mesh_trial_ref,
					 double* mesh_grad_trial_ref,
					 double* mesh_dof,
					 int* mesh_l2g,
					 double* dV_ref,
					 double* u_trial_ref,
					 double* u_grad_trial_ref,
					 double* u_test_ref,
					 double* u_grad_test_ref,
					 //element boundary
					 double* mesh_trial_trace_ref,
					 double* mesh_grad_trial_trace_ref,
					 double* dS_ref,
					 double* u_trial_trace_ref,
					 double* u_grad_trial_trace_ref,
					 double* u_test_trace_ref,
					 double* u_grad_test_trace_ref,
					 double* normal_ref,
					 double* boundaryJac_ref,
					 //physics
					 int nElements_global,
					 int* u_l2g,
					 double* u_dof,
					 int* sd_rowptr,
					 int* sd_colind,
					 double* q_a,
					 double* q_r,
					 double* elementJacobian_u_u,
					 double* element_u,
					 int eN)
    {
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    elementJacobian_u_u[i*nDOF_trial_element+j]=0.0;
	  }
      for  (int k=0;k<nQuadraturePoints_element;k++)
	{
	  int eN_k = eN*nQuadraturePoints_element+k; //index to a scalar at a quadrature point
	  
	  //declare local storage
	  register double u=0.0,
	    grad_u[nSpace],
	    *a=NULL,
	    dr=0.0,
	    jac[nSpace*nSpace],
	    jacDet,
	    jacInv[nSpace*nSpace],
	    u_grad_trial[nDOF_trial_element*nSpace],
	    dV,
	    u_test_dV[nDOF_test_element],
	    u_grad_test_dV[nDOF_test_element*nSpace],
	    x,y,z,
	    G[nSpace*nSpace],G_dd_G,tr_G;
	  //
	  //calculate solution and gradients at quadrature points
	  //
	  ck.calculateMapping_element(eN,
				      k,
				      mesh_dof,
				      mesh_l2g,
				      mesh_trial_ref,
				      mesh_grad_trial_ref,
				      jac,
				      jacDet,
				      jacInv,
				      x,y,z);
	  //get the physical integration weight
	  dV = fabs(jacDet)*dV_ref[k];
	  //get metric tensor and friends
	  ck.calculateG(jacInv,G,G_dd_G,tr_G);
	  //get the trial function gradients
	  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	  //get the solution 	
	  ck.valFromElementDOF(element_u,&u_trial_ref[k*nDOF_trial_element],u);
	  //get the solution gradients
	  ck.gradFromElementDOF(element_u,u_grad_trial,grad_u);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
	      for (int I=0;I<nSpace;I++)
		{
		  u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	    }
	  //
	  //calculate pde coefficients and derivatives at quadrature points
	  //
	  //evaluateCoefficients()
	  a = &q_a[eN_k*sd_rowptr[nSpace]];
	  dr = 0.0;
	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int i_nSpace=i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int j_nSpace = j*nSpace;
		  elementJacobian_u_u[i*nDOF_trial_element+j] += ck.SimpleDiffusionJacobian_weak(sd_rowptr,sd_colind,a,&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]) +
		    ck.ReactionJacobian_weak(dr,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]); 		     
		}//j
	    }//i
	}//k
    }
    void calculateJacobian(//element
			   double* mesh_trial_ref,
			   double* mesh_grad_trial_ref,
			   double* mesh_dof,
			   int* mesh_l2g,
			   double* dV_ref,
			   double* u_trial_ref,
			   double* u_grad_trial_ref,
			   double* u_test_ref,
			   double* u_grad_test_ref,
			   //element boundary
			   double* mesh_trial_trace_ref,
			   double* mesh_grad_trial_trace_ref,
			   double* dS_ref,
			   double* u_trial_trace_ref,
			   double* u_grad_trial_trace_ref,
			   double* u_test_trace_ref,
			   double* u_grad_test_trace_ref,
			   double* normal_ref,
			   double* boundaryJac_ref,
			   //physics
			   int nElements_global,
			   int* u_l2g,
			   double* u_dof,
			   int* sd_rowptr,
			   int* sd_colind,
			   double* q_a,
			   double* q_r,
			   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			   double* globalJacobian,
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
			   double* ebqe_a,
			   int* isDOFBoundary_u,
			   double* ebqe_bc_u_ext,
			   int* isDiffusiveFluxBoundary_u,
			   double* ebqe_bc_flux_u_ext,
			   int* csrColumnOffsets_eb_u_u,
			   double* ebqe_penalty_ext)
    {
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element*nDOF_trial_element],element_u[nDOF_trial_element];
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      register int eN_j = eN*nDOF_trial_element+j;
	      element_u[j] = u_dof[u_l2g[eN_j]];
	    }
	  calculateElementJacobian(mesh_trial_ref,
				   mesh_grad_trial_ref,
				   mesh_dof,
				   mesh_l2g,
				   dV_ref,
				   u_trial_ref,
				   u_grad_trial_ref,
				   u_test_ref,
				   u_grad_test_ref,
				   mesh_trial_trace_ref,
				   mesh_grad_trial_trace_ref,
				   dS_ref,
				   u_trial_trace_ref,
				   u_grad_trial_trace_ref,
				   u_test_trace_ref,
				   u_grad_test_trace_ref,
				   normal_ref,
				   boundaryJac_ref,
				   nElements_global,
				   u_l2g,
				   u_dof,
				   sd_rowptr,
				   sd_colind,
				   q_a,
				   q_r,
				   elementJacobian_u_u,
				   element_u,
				   eN);
	  //
	  //load into element Jacobian into global Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  int eN_i_j = eN_i*nDOF_trial_element+j;
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i*nDOF_trial_element+j];
		}//j
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
      //
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0],
            ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
            eN_nDOF_trial_element = eN*nDOF_trial_element;
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;

	      register double u_ext=0.0,
		grad_u_ext[nSpace],
		m_ext=0.0,
		dm_ext=0.0,
		*a_ext,
		//		f_ext[nSpace],
		//df_ext[nSpace],
		r_ext=0.0,
		dflux_u_u_ext=0.0,
		bc_u_ext=0.0,
		//bc_grad_u_ext[nSpace],
		bc_m_ext=0.0,
		bc_dm_ext=0.0,
		bc_f_ext[nSpace],
		bc_df_ext[nSpace],
		fluxJacobian_u_u[nDOF_trial_element],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		u_grad_trial_trace[nDOF_trial_element*nSpace],
		normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		//
		G[nSpace*nSpace],G_dd_G,tr_G;
	      // 
	      //calculate the solution and gradients at quadrature points 
	      // 
	      // u_ext=0.0;
	      // for (int I=0;I<nSpace;I++)
	      //   {
	      //     grad_u_ext[I] = 0.0;
	      //     bc_grad_u_ext[I] = 0.0;
	      //   }
	      // for (int j=0;j<nDOF_trial_element;j++) 
	      //   { 
	      //     register int eN_j = eN*nDOF_trial_element+j,
	      //       ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
	      //       ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	      //     u_ext += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial_ext[ebNE_kb_j]); 
	                     
	      //     for (int I=0;I<nSpace;I++)
	      //       {
	      //         grad_u_ext[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
	      //       } 
	      //   }
	      ck.calculateMapping_elementBoundary(eN,
						  ebN_local,
						  kb,
						  ebN_local_kb,
						  mesh_dof,
						  mesh_l2g,
						  mesh_trial_trace_ref,
						  mesh_grad_trial_trace_ref,
						  boundaryJac_ref,
						  jac_ext,
						  jacDet_ext,
						  jacInv_ext,
						  boundaryJac,
						  metricTensor,
						  metricTensorDetSqrt,
						  normal_ref,
						  normal,
						  x_ext,y_ext,z_ext);
	      /* ck.calculateMappingVelocity_elementBoundary(eN, */
	      /* 						  ebN_local, */
	      /* 						  kb, */
	      /* 						  ebN_local_kb, */
	      /* 						  mesh_velocity_dof, */
	      /* 						  mesh_l2g, */
	      /* 						  mesh_trial_trace_ref, */
	      /* 						  xt_ext,yt_ext,zt_ext, */
	      /* 						  normal, */
	      /* 						  boundaryJac, */
	      /* 						  metricTensor, */
	      /* 						  integralScaling); */
	      /* dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb]; */
	      dS = metricTensorDetSqrt*dS_ref[kb];
	      //dS = metricTensorDetSqrt*dS_ref[kb];
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial_trace,grad_u_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		}
	      //
	      //load the boundary values
	      //
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	      // 
	      //calculate the internal and external trace of the pde coefficients 
	      // 
	      /* evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace], */
	      /* 			   u_ext, */
	      /* 			   //VRANS */
	      /* 			   porosity_ext, */
	      /* 			   // */
	      /* 			   m_ext, */
	      /* 			   dm_ext, */
	      /* 			   f_ext, */
	      /* 			   df_ext); */
	      /* evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace], */
	      /* 			   bc_u_ext, */
	      /* 			   //VRANS */
	      /* 			   porosity_ext, */
	      /* 			   // */
	      /* 			   bc_m_ext, */
	      /* 			   bc_dm_ext, */
	      /* 			   bc_f_ext, */
	      /* 			   bc_df_ext); */
	      a_ext = &ebqe_a[ebNE_kb*sd_rowptr[nSpace]];
	      //
	      //moving domain
	      //
	      /* double velocity_ext[nSpace]; */
	      /* double mesh_velocity[3]; */
	      /* for (int I=0;I<nSpace;I++) */
	      /* 	velocity_ext[I] = ebqe_velocity_ext[ebNE_kb_nSpace+0];// - MOVING_DOMAIN*mesh_velocity[I]; */
	      // 
	      //calculate the numerical fluxes 
	      // 
	      /* exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u[ebNE_kb], */
	      /* 					       isFluxBoundary_u[ebNE_kb], */
	      /* 					       normal, */
	      /* 					       df_ext, */
	      /* 					       dflux_u_u_ext); */
	      //
	      //calculate the flux jacobian
	      //
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
		  register int j_nSpace = j*nSpace,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
		  fluxJacobian_u_u[j]=ExteriorNumericalDiffusiveFluxJacobian(sd_rowptr,
									     sd_colind,
									     isDOFBoundary_u[ebNE_kb],
									     isDiffusiveFluxBoundary_u[ebNE_kb],
									     normal,
									     a_ext,
									     u_trial_trace_ref[ebN_local_kb_j],
									     &u_grad_trial_trace[j_nSpace],
									     ebqe_penalty_ext[ebNE_kb]);
		  //ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref[ebN_local_kb_j]);
		}//j
	      //
	      //update the global Jacobian from the flux Jacobian
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i = eN*nDOF_test_element+i;
		  //register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;

		      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*u_test_dS[i];
		    }//j
		}//i
	    }//kb
	}//ebNE
    }//computeJacobian
  };//cppADR

  inline cppADR_base* newADR(int nSpaceIn,
			     int nQuadraturePoints_elementIn,
			     int nDOF_mesh_trial_elementIn,
			     int nDOF_trial_elementIn,
			     int nDOF_test_elementIn,
			     int nQuadraturePoints_elementBoundaryIn,
			     int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<cppADR_base,cppADR,CompKernel>(nSpaceIn,
										       nQuadraturePoints_elementIn,
										       nDOF_mesh_trial_elementIn,
										       nDOF_trial_elementIn,
										       nDOF_test_elementIn,
										       nQuadraturePoints_elementBoundaryIn,
										       CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<cppADR_base,cppADR,CompKernel>(nSpaceIn,
										     nQuadraturePoints_elementIn,
										     nDOF_mesh_trial_elementIn,
										     nDOF_trial_elementIn,
										     nDOF_test_elementIn,
										     nQuadraturePoints_elementBoundaryIn,
										     CompKernelFlag);
  }
}//proteus
#endif
