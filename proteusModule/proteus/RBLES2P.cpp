#include "RBLES2P.h"
#include <iostream>

extern "C" void calculateResidual_RBLES2P(int nElements_global,
					 double alphaBDF,
					 double epsFact_rho,
					 double epsFact_mu, 
					 double sigma,
					 double rho_0,
					 double nu_0,
					 double rho_1,
					 double nu_1,
					 double hFactor,
					 double shockCapturingDiffusion,
					 int* p_l2g, int* vel_l2g, 
					 double* elementDiameter,
					 double* p_dof, double* u_dof, double* v_dof, double* w_dof,
					 double* p_trial, double* vel_trial,
					 double* p_grad_trial, double* vel_grad_trial,
					 double* p_test_dV, double* vel_test_dV,
					 double* p_grad_test_dV, double* vel_grad_test_dV,
					 double* vel_Hess_trial,double* vel_Hess_test_dV,
					 double* g,
					 double* phi,
					 double* n,
					 double* kappa,
					 double* q_mom_u_acc,
					 double* q_mom_v_acc,
					 double* q_mom_w_acc,
					 double* q_mass_adv,
					 double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
					 double* q_velocity_last,
					 double* q_cfl,
					 double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
					 double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
					 double* q_elementResidual_p, double* q_elementResidual_u, double* q_elementResidual_v, double* q_elementResidual_w,
					 int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
					 int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
					 int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
					 int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
					 int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
					 int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
					 int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
					 int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
					 int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
					 int offset_p, int offset_u, int offset_v, int offset_w, int stride_p, 
					 int stride_u, int stride_v, int stride_w, double* globalResidual,
					 int nExteriorElementBoundaries_global,
					 int* exteriorElementBoundariesArray,
					 int* elementBoundaryElementsArray,
					 int* elementBoundaryLocalElementBoundariesArray,
					 double* p_trial_ext,
					 double* vel_trial_ext,
					 double* p_grad_trial_ext, 
					 double* vel_grad_trial_ext,
					 double* ebqe_phi_ext,
					 double* ebqe_n_ext,
					 double* ebqe_kappa_ext,
					 int* isDOFBoundary_p,
					 int* isDOFBoundary_u,
					 int* isDOFBoundary_v,
					 int* isDOFBoundary_w,
					 int* isAdvectiveFluxBoundary_p,
					 int* isAdvectiveFluxBoundary_u,
					 int* isAdvectiveFluxBoundary_v,
					 int* isAdvectiveFluxBoundary_w,
					 int* isDiffusiveFluxBoundary_u,
					 int* isDiffusiveFluxBoundary_v,
					 int* isDiffusiveFluxBoundary_w,
					 double* ebqe_bc_p_ext,
					 double* ebqe_bc_flux_mass_ext,
					 double* ebqe_bc_flux_mom_u_adv_ext,
					 double* ebqe_bc_flux_mom_v_adv_ext,
					 double* ebqe_bc_flux_mom_w_adv_ext,
					 double* ebqe_bc_u_ext,
					 double* ebqe_bc_flux_u_diff_ext,
					 double* ebqe_penalty_ext,
					 double* ebqe_bc_v_ext,
					 double* ebqe_bc_flux_v_diff_ext,
					 double* ebqe_bc_w_ext,
					 double* ebqe_bc_flux_w_diff_ext,
					 double* p_test_dS_ext,
					 double* vel_test_dS_ext,
					 double* q_velocity,
					 double* ebqe_velocity,
					 double* flux)
{
  using namespace RBLES2P;
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
  double globalConservationError=0.0;
  for(int eN=0;eN<nElements_global;eN++)
    {
      //declare local storage for element residual and initialize
      register double elementResidual_p[nDOF_test_element],
	elementResidual_u[nDOF_test_element],
	elementResidual_v[nDOF_test_element],
	elementResidual_w[nDOF_test_element];
      const double eps_rho = epsFact_rho*elementDiameter[eN],
      	eps_mu = epsFact_mu*elementDiameter[eN];
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_p[i]=0.0;
	  elementResidual_u[i]=0.0;
	  elementResidual_v[i]=0.0;
	  elementResidual_w[i]=0.0;
	}//i
      //loop over quadrature points and compute integrands
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  //compute indeces and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k,
	    eN_k_nSpace = eN_k*nSpace;
	  register double p=0.0,u=0.0,v=0.0,w=0.0,grad_p[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
	    mom_u_acc=0.0,
	    dmom_u_acc_u=0.0,
	    mom_v_acc=0.0,
	    dmom_v_acc_v=0.0,
	    mom_w_acc=0.0,
	    dmom_w_acc_w=0.0,
	    mom_u_acc_t=0.0,
	    dmom_u_acc_u_t=0.0,
	    mom_v_acc_t=0.0,
	    dmom_v_acc_v_t=0.0,
	    mom_w_acc_t=0.0,
	    dmom_w_acc_w_t=0.0,
	    pdeResidual_p=0.0,
	    pdeResidual_u=0.0,
	    pdeResidual_v=0.0,
	    pdeResidual_w=0.0,
	    subgrid_p=0.0,
	    subgrid_u=0.0,
	    subgrid_v=0.0,
	    subgrid_w=0.0,
	    tau_0=0.0,
	    tau_1=0.0;
	    
	 // double& N = vel_trial[eN_k_j]  

 

          /*for (int j=0;j<nDOF_trial_element;j++)
            {
              int eN_k_j=eN_k*nDOF_trial_element+j;
	      Nq[j] =   p_test_dV [eN_k_j];
	      Nw[j] = vel_test_dV [eN_k_j];
              Np[j] =   p_trial[eN_k_j];
	      Nu[j] = vel_trial[eN_k_j];

              eN_k_j*=nSpace;
	    for (int I=0;I<nSpace;I++)
	      {                
		gNq[j][I] =   p_grad_test_dV [eN_k_j];
	        gNw[j][I] = vel_grad_test_dV [eN_k_j];
                gNp[j][I] =   p_grad_trial[eN_k_j];
	        gNu[j][I] = vel_grad_trial[eN_k_j];
                eN_k_j += 1;
	      }
            }*/


	    
          //
          //compute solution and gradients at quadrature points
          //
	  p=0.0;u=0.0;v=0.0;w=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_p[I]=0.0;grad_u[I]=0.0;grad_v[I]=0.0;grad_w[I]=0.0;
	    }
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_j=eN*nDOF_trial_element+j;
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
              p += valFromDOF_c(p_dof[  p_l2g[eN_j]],  p_trial[eN_k_j]);
              u += valFromDOF_c(u_dof[vel_l2g[eN_j]],vel_trial[eN_k_j]);
              v += valFromDOF_c(v_dof[vel_l2g[eN_j]],vel_trial[eN_k_j]);
              w += valFromDOF_c(w_dof[vel_l2g[eN_j]],vel_trial[eN_k_j]);
	      
	      for (int I=0;I<nSpace;I++)
		{
		  grad_p[I] += gradFromDOF_c(p_dof[  p_l2g[eN_j]],  p_grad_trial[eN_k_j_nSpace+I]);
		  grad_u[I] += gradFromDOF_c(u_dof[vel_l2g[eN_j]],vel_grad_trial[eN_k_j_nSpace+I]);
		  grad_v[I] += gradFromDOF_c(v_dof[vel_l2g[eN_j]],vel_grad_trial[eN_k_j_nSpace+I]);
		  grad_w[I] += gradFromDOF_c(w_dof[vel_l2g[eN_j]],vel_grad_trial[eN_k_j_nSpace+I]);
		}
	    }
	  //save velocity at quadrature points for other models to use
	  // q_velocity[eN_k_nSpace+0]=u+subgridError_u;
	  // q_velocity[eN_k_nSpace+1]=v+subgridError_v;
	  // q_velocity[eN_k_nSpace+2]=w+subgridError_w;
// 	  q_velocity[eN_k_nSpace+0]=u;
// 	  q_velocity[eN_k_nSpace+1]=v;
// 	  q_velocity[eN_k_nSpace+2]=w;
          //
          //calculate pde coefficients at quadrature points
          //
 
          double rho,nu,mu,H_rho,d_rho,H_mu,d_mu,norm_n;
          H_rho = smoothedHeaviside(eps_rho,phi[eN_k]);
          d_rho = smoothedDirac(eps_rho,phi[eN_k]);
          H_mu = smoothedHeaviside(eps_mu,phi[eN_k]);
          d_mu = smoothedDirac(eps_mu,phi[eN_k]);
  
          rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
          nu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
          mu  = rho_0*nu_0*(1.0-H_mu)+rho_1*nu_1*H_mu;
 
           //u momentum accumulation
          mom_u_acc=u;
          dmom_u_acc_u=1.0;
  
          //v momentum accumulation
          mom_v_acc=v;
          dmom_v_acc_v=1.0;
  
          //w momentum accumulation
          mom_w_acc=w;
          dmom_w_acc_w=1.0;
          
	  //
	  //save momentum for time history 
	  //
	  q_mom_u_acc[eN_k] = mom_u_acc;			    
	  q_mom_v_acc[eN_k] = mom_v_acc;			    
	  q_mom_w_acc[eN_k] = mom_w_acc;
          
          //
          //calculate time derivative at quadrature points
          //
          bdf_c(alphaBDF,
		q_mom_u_acc_beta_bdf[eN_k],
		mom_u_acc,
		dmom_u_acc_u,
		mom_u_acc_t,
		dmom_u_acc_u_t);
          bdf_c(alphaBDF,
		q_mom_v_acc_beta_bdf[eN_k],
		mom_v_acc,
		dmom_v_acc_v,
		mom_v_acc_t,
		dmom_v_acc_v_t);
          bdf_c(alphaBDF,
		q_mom_w_acc_beta_bdf[eN_k],
		mom_w_acc,
		dmom_w_acc_w,
		mom_w_acc_t,
		dmom_w_acc_w_t);          
          
         
          //
          //calculate subgrid error - strong residual
          //
          //calculate strong residual
	  pdeResidual_p = grad_u[0] + grad_v[1] + grad_w[2];

	  pdeResidual_u = rho *(mom_u_acc_t + u*grad_u[0] + v*grad_u[1] + w*grad_u[2] - g[0])
	                + grad_p[0]; // - mu poisson_u  ==> ommited as only first order right now !! 

	  pdeResidual_v = rho *(mom_v_acc_t + u*grad_v[0] + v*grad_v[1] + w*grad_v[2] - g[1])
	                + grad_p[1]; // - mu poisson_u  ==> ommited as only first order right now !! 
	                
	  pdeResidual_w = rho *(mom_w_acc_t + u*grad_w[0] + v*grad_w[1] + w*grad_w[2] - g[2])
	                + grad_p[2]; // - mu poisson_u  ==> ommited as only first order right now !! 

          //calculate tau and subgrid 
          calculateSubgridError_tau_c(hFactor,elementDiameter[eN],
				      dmom_u_acc_u_t,rho,mu,
				      u,v,w,
				      tau_0,tau_1,q_cfl[eN_k]);



          subgrid_u = -tau_0*pdeResidual_u;
          subgrid_v = -tau_0*pdeResidual_v;
          subgrid_w = -tau_0*pdeResidual_w;

          subgrid_p = -tau_1*pdeResidual_p;

          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_k_i=eN_k*nDOF_test_element+i,
		eN_k_i_nSpace = eN_k_i*nSpace;

	      elementResidual_p[i] += pdeResidual_p*p_test_dV[eN_k_i]
	                            - subgrid_u*p_grad_test_dV[eN_k_i_nSpace+0]
	                            - subgrid_v*p_grad_test_dV[eN_k_i_nSpace+1]
	                            - subgrid_w*p_grad_test_dV[eN_k_i_nSpace+2];

	      elementResidual_u[i] += rho *(mom_u_acc_t + u*grad_u[0] + v*grad_u[1] + w*grad_u[2] - g[0])*vel_test_dV[eN_k_i]
	                            + mu*(grad_u[0]*vel_grad_test_dV[eN_k_i_nSpace+0]
	                                 +grad_u[1]*vel_grad_test_dV[eN_k_i_nSpace+1]
	                                 +grad_u[2]*vel_grad_test_dV[eN_k_i_nSpace+2])	
	                            -p*vel_grad_test_dV[eN_k_i_nSpace+0]     	 
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+0] * (subgrid_u*u)// + u*subgrid_u + subgrid_u*subgrid_u)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+1] * (subgrid_u*v)// + u*subgrid_v + subgrid_u*subgrid_v)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+2] * (subgrid_u*w);// + u*subgrid_w + subgrid_u*subgrid_w);
		            
	      elementResidual_v[i] += rho *(mom_v_acc_t + u*grad_v[0] + v*grad_v[1] + w*grad_v[2] - g[1])*vel_test_dV[eN_k_i]
	                            + mu*(grad_v[0]*vel_grad_test_dV[eN_k_i_nSpace+0]
	                                 +grad_v[1]*vel_grad_test_dV[eN_k_i_nSpace+1]
	                                 +grad_v[2]*vel_grad_test_dV[eN_k_i_nSpace+2])		
	                            -p*vel_grad_test_dV[eN_k_i_nSpace+1] 	                                  
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+0] * (subgrid_v*u)// + v*subgrid_u + subgrid_v*subgrid_u)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+1] * (subgrid_v*v)// + v*subgrid_v + subgrid_v*subgrid_v)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+2] * (subgrid_v*w);// + v*subgrid_w + subgrid_v*subgrid_w);


	      elementResidual_w[i] += rho *(mom_w_acc_t + u*grad_w[0] + v*grad_w[1] + w*grad_w[2] - g[2])*vel_test_dV[eN_k_i]
	                            + mu*(grad_w[0]*vel_grad_test_dV[eN_k_i_nSpace+0]
	                                 +grad_w[1]*vel_grad_test_dV[eN_k_i_nSpace+1]
	                                 +grad_w[2]*vel_grad_test_dV[eN_k_i_nSpace+2])		 
	                            -p*vel_grad_test_dV[eN_k_i_nSpace+2] 
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+0] * (subgrid_w*u)// + w*subgrid_u + subgrid_w*subgrid_u)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+1] * (subgrid_w*v)// + w*subgrid_v + subgrid_w*subgrid_v)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+2] * (subgrid_w*w);// + w*subgrid_w + subgrid_w*subgrid_w);
            }//i
	}
      //
      //load element into global residual and save element residual
      //
      for(int i=0;i<nDOF_test_element;i++) 
        { 
          register int eN_i=eN*nDOF_test_element+i;
          
          q_elementResidual_p[eN_i]+=elementResidual_p[i];
          q_elementResidual_u[eN_i]+=elementResidual_u[i];
          q_elementResidual_v[eN_i]+=elementResidual_v[i];
          q_elementResidual_w[eN_i]+=elementResidual_w[i];

          globalResidual[offset_p+stride_p*  p_l2g[eN_i]]+=elementResidual_p[i];
          globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
          globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
          globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
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
	eN  = elementBoundaryElementsArray[ebN*2+0];
      register double elementResidual_p[nDOF_test_element],
	              elementResidual_u[nDOF_test_element],
	              elementResidual_v[nDOF_test_element],
	              elementResidual_w[nDOF_test_element];
      const double eps_rho = epsFact_rho*elementDiameter[eN],
      	           eps_mu  = epsFact_mu *elementDiameter[eN];
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_p[i]=0.0;
	  elementResidual_u[i]=0.0;
	  elementResidual_v[i]=0.0;
	  elementResidual_w[i]=0.0;
	}
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  register double p_ext=0.0,
	    u_ext=0.0,
	    v_ext=0.0,
	    w_ext=0.0,
	    gamma=0.0,	    	    
	    grad_p_ext[nSpace],
	    grad_u_ext[nSpace],
	    grad_v_ext[nSpace],
	    grad_w_ext[nSpace];
           double norm[nSpace];
	   	
	  norm[nSpace] = ebqe_n_ext[ebNE_kb_nSpace];    
	  // 
	  //calculate the solution and gradients at quadrature points 
	  // 
	  p_ext=0.0;u_ext=0.0;v_ext=0.0;w_ext=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_p_ext[I] = 0.0;
	      grad_u_ext[I] = 0.0;
	      grad_v_ext[I] = 0.0;
	      grad_w_ext[I] = 0.0;
	    }
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      int eN_j = eN*nDOF_trial_element+j;
	      int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
	      int ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	      p_ext += valFromDOF_c(p_dof[p_l2g[eN_j]],p_trial_ext[ebNE_kb_j]); 
	      u_ext += valFromDOF_c(u_dof[vel_l2g[eN_j]],vel_trial_ext[ebNE_kb_j]); 
	      v_ext += valFromDOF_c(v_dof[vel_l2g[eN_j]],vel_trial_ext[ebNE_kb_j]); 
	      w_ext += valFromDOF_c(w_dof[vel_l2g[eN_j]],vel_trial_ext[ebNE_kb_j]); 
               
	      for (int I=0;I<nSpace;I++)
		{
		  grad_p_ext[I] += gradFromDOF_c(p_dof[p_l2g[eN_j]],p_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_u_ext[I] += gradFromDOF_c(u_dof[vel_l2g[eN_j]],vel_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_v_ext[I] += gradFromDOF_c(v_dof[vel_l2g[eN_j]],vel_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_w_ext[I] += gradFromDOF_c(w_dof[vel_l2g[eN_j]],vel_grad_trial_ext[ebNE_kb_j_nSpace+I]);
		} 
	    }
	    
	    
          double rho,nu,mu,H_rho,d_rho,H_mu,d_mu,norm_n;
          H_rho = smoothedHeaviside(eps_rho,phi[eN]);
          d_rho = smoothedDirac(eps_rho,phi[eN]);
          H_mu = smoothedHeaviside(eps_mu,phi[eN]);
          d_mu = smoothedDirac(eps_mu,phi[eN]);
  
          rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
          nu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
          mu  = rho_0*nu_0*(1.0-H_mu)+rho_1*nu_1*H_mu;
	            
          calculateInteriorPenalty   (hFactor,elementDiameter[eN],
				      mu,gamma);


          for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
	      int ebNE_kb_i_nSpace= ebNE_kb_i*nSpace;
	      elementResidual_p[i]+= (u_ext*norm[0] + v_ext*norm[1] + w_ext*norm[2])*p_test_dS_ext[ebNE_kb_i];
	      elementResidual_u[i]+= (p_ext*norm[0] - mu*(grad_u_ext[0]*norm[0] 
	                                                + grad_u_ext[1]*norm[1] 
							+ grad_u_ext[2]*norm[2]))*vel_test_dS_ext[ebNE_kb_i]							
				// - mu*(u_ext)*(vel_grad_test_dS_ext[ebNE_kb_i_nSpace+0]*norm[0] 
  				//	       + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+1]*norm[1] 
  				//	       + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+2]*norm[2])
			           + gamma*(u_ext)*vel_test_dS_ext[ebNE_kb_i];		  
							
	      elementResidual_v[i]+= (p_ext*norm[1] - mu*(grad_v_ext[0]*norm[0] 
	                                                + grad_v_ext[1]*norm[1] 
							+ grad_v_ext[2]*norm[2]))*vel_test_dS_ext[ebNE_kb_i]
				// - mu*(v_ext)*(vel_grad_test_dS_ext[ebNE_kb_i_nSpace+0]*norm[0] 
  				//	       + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+1]*norm[1] 
  				//	       + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+2]*norm[2])
                                   + gamma*(v_ext)*vel_test_dS_ext[ebNE_kb_i];	
				   				       
	      elementResidual_w[i]+= (p_ext*norm[2] - mu*(grad_w_ext[0]*norm[0] 
	                                                + grad_w_ext[1]*norm[1] 
							+ grad_w_ext[2]*norm[2]))*vel_test_dS_ext[ebNE_kb_i]
				 //- mu*(w_ext)*(vel_grad_test_dS_ext[ebNE_kb_i_nSpace+0]*norm[0] 
  				 //	       + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+1]*norm[1] 
  				 //	       + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+2]*norm[2])
                                   + gamma*(w_ext)*vel_test_dS_ext[ebNE_kb_i];
	    }//i
	    


 
	}//kb
      //
      //update the element and global residual storage
      //
      for (int i=0;i<nDOF_test_element;i++)
	{
	  int eN_i = eN*nDOF_test_element+i;
	  q_elementResidual_p[eN_i]+=elementResidual_p[i];
	  q_elementResidual_u[eN_i]+=elementResidual_u[i];
	  q_elementResidual_v[eN_i]+=elementResidual_v[i];
	  q_elementResidual_w[eN_i]+=elementResidual_w[i];

	  globalResidual[offset_p+stride_p*  p_l2g[eN_i]]+=elementResidual_p[i];
	  globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
	  globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
	  globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
	}//i
    }//ebNE


}

extern "C" void calculateJacobian_RBLES2P(int nElements_global,
					 double alphaBDF,
					 double epsFact_rho,
					 double epsFact_mu,
					 double sigma,
					 double rho_0,
					 double nu_0,
					 double rho_1,
					 double nu_1,
					 double hFactor,
					 double shockCapturingDiffusion,
					 int* p_l2g, int* vel_l2g,
					 double* elementDiameter,
					 double* p_dof, double* u_dof, double* v_dof, double* w_dof,
					 double* p_trial, double* vel_trial,
					 double* p_grad_trial, double* vel_grad_trial,
					 double* p_test_dV, double* vel_test_dV, 
					 double* p_grad_test_dV, double* vel_grad_test_dV,
					 double* vel_Hess_trial,double* vel_Hess_test_dV,
					 double* g,
					 double* phi,
					 double* n,
					 double* kappa,
					 double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
					 double* q_velocity_last,
					 double* q_cfl,
					 double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
					 int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
					 int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
					 int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
					 int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
					 int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
					 int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
					 int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
					 int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
					 int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
					 int* csrRowIndeces_p_p,int* csrColumnOffsets_p_p,
					 int* csrRowIndeces_p_u,int* csrColumnOffsets_p_u,
					 int* csrRowIndeces_p_v,int* csrColumnOffsets_p_v,
					 int* csrRowIndeces_p_w,int* csrColumnOffsets_p_w,
					 int* csrRowIndeces_u_p,int* csrColumnOffsets_u_p,
					 int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
					 int* csrRowIndeces_u_v,int* csrColumnOffsets_u_v,
					 int* csrRowIndeces_u_w,int* csrColumnOffsets_u_w,
					 int* csrRowIndeces_v_p,int* csrColumnOffsets_v_p,
					 int* csrRowIndeces_v_u,int* csrColumnOffsets_v_u,
					 int* csrRowIndeces_v_v,int* csrColumnOffsets_v_v,
					 int* csrRowIndeces_v_w,int* csrColumnOffsets_v_w,
					 int* csrRowIndeces_w_p,int* csrColumnOffsets_w_p,
					 int* csrRowIndeces_w_u,int* csrColumnOffsets_w_u,
					 int* csrRowIndeces_w_v,int* csrColumnOffsets_w_v,
					 int* csrRowIndeces_w_w,int* csrColumnOffsets_w_w,
					 double* globalJacobian,
					 int nExteriorElementBoundaries_global,
					 int* exteriorElementBoundariesArray,
					 int* elementBoundaryElementsArray,
					 int* elementBoundaryLocalElementBoundariesArray,
					 double* p_trial_ext,
					 double* vel_trial_ext,
					 double* p_grad_trial_ext,
					 double* vel_grad_trial_ext,
					 double* ebqe_phi_ext,
					 double* ebqe_n_ext,
					 double* ebqe_kappa_ext,
					 int* isDOFBoundary_p,
					 int* isDOFBoundary_u,
					 int* isDOFBoundary_v,
					 int* isDOFBoundary_w,
					 int* isAdvectiveFluxBoundary_p,
					 int* isAdvectiveFluxBoundary_u,
					 int* isAdvectiveFluxBoundary_v,
					 int* isAdvectiveFluxBoundary_w,
					 int* isDiffusiveFluxBoundary_u,
					 int* isDiffusiveFluxBoundary_v,
					 int* isDiffusiveFluxBoundary_w,
					 double* ebqe_bc_p_ext,
					 double* ebqe_bc_flux_mass_ext,
					 double* ebqe_bc_flux_mom_u_adv_ext,
					 double* ebqe_bc_flux_mom_v_adv_ext,
					 double* ebqe_bc_flux_mom_w_adv_ext,
					 double* ebqe_bc_u_ext,
					 double* ebqe_bc_flux_u_diff_ext,
					 double* ebqe_penalty_ext,
					 double* ebqe_bc_v_ext,
					 double* ebqe_bc_flux_v_diff_ext,
					 double* ebqe_bc_w_ext,
					 double* ebqe_bc_flux_w_diff_ext,
					 double* p_test_dS_ext,
					 double* vel_test_dS_ext,
					 int* csrColumnOffsets_eb_p_p,
					 int* csrColumnOffsets_eb_p_u,
					 int* csrColumnOffsets_eb_p_v,
					 int* csrColumnOffsets_eb_p_w,
					 int* csrColumnOffsets_eb_u_p,
					 int* csrColumnOffsets_eb_u_u,
					 int* csrColumnOffsets_eb_u_v,
					 int* csrColumnOffsets_eb_u_w,
					 int* csrColumnOffsets_eb_v_p,
					 int* csrColumnOffsets_eb_v_u,
					 int* csrColumnOffsets_eb_v_v,
					 int* csrColumnOffsets_eb_v_w,
					 int* csrColumnOffsets_eb_w_p,
					 int* csrColumnOffsets_eb_w_u,
					 int* csrColumnOffsets_eb_w_v,
					 int* csrColumnOffsets_eb_w_w)
{
  using namespace RBLES2P;
  //
  //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
  //
  for(int eN=0;eN<nElements_global;eN++)
    {
      const double eps_rho = epsFact_rho*elementDiameter[eN],
      	eps_mu = epsFact_mu*elementDiameter[eN];

      register double  elementJacobian_p_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_p_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_p_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_p_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_p[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_w[nDOF_test_element][nDOF_trial_element];
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    elementJacobian_p_p[i][j]=0.0;
	    elementJacobian_p_u[i][j]=0.0;
	    elementJacobian_p_v[i][j]=0.0;
	    elementJacobian_p_w[i][j]=0.0;
	    elementJacobian_u_p[i][j]=0.0;
	    elementJacobian_u_u[i][j]=0.0;
	    elementJacobian_u_v[i][j]=0.0;
	    elementJacobian_u_w[i][j]=0.0;
	    elementJacobian_v_p[i][j]=0.0;
	    elementJacobian_v_u[i][j]=0.0;
	    elementJacobian_v_v[i][j]=0.0;
	    elementJacobian_v_w[i][j]=0.0;
	    elementJacobian_w_p[i][j]=0.0;
	    elementJacobian_w_u[i][j]=0.0;
	    elementJacobian_w_v[i][j]=0.0;
	    elementJacobian_w_w[i][j]=0.0;
	  }
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
	    eN_k_nSpace = eN_k*nSpace; //index to a vector at a quadrature point

	  //declare local storage
	  register double p=0.0,u=0.0,v=0.0,w=0.0,
	    grad_p[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
	    mom_u_acc=0.0,
	    dmom_u_acc_u=0.0,
	    mom_v_acc=0.0,
	    dmom_v_acc_v=0.0,
	    mom_w_acc=0.0,
	    dmom_w_acc_w=0.0,
	    mom_u_acc_t=0.0,
	    dmom_u_acc_u_t=0.0,
	    mom_v_acc_t=0.0,
	    dmom_v_acc_v_t=0.0,
	    mom_w_acc_t=0.0,
	    dmom_w_acc_w_t=0.0,
	    pdeResidual_p=0.0,
	    pdeResidual_u=0.0,
	    pdeResidual_v=0.0,
	    pdeResidual_w=0.0,
	    subgrid_p=0.0,
	    subgrid_u=0.0,
	    subgrid_v=0.0,
	    subgrid_w=0.0,
	    dpdeResidual_p_u[nDOF_trial_element],dpdeResidual_p_v[nDOF_trial_element],dpdeResidual_p_w[nDOF_trial_element],
	    dpdeResidual_u_p[nDOF_trial_element],dpdeResidual_u_u[nDOF_trial_element],
	    dpdeResidual_v_p[nDOF_trial_element],dpdeResidual_v_v[nDOF_trial_element],
	    dpdeResidual_w_p[nDOF_trial_element],dpdeResidual_w_w[nDOF_trial_element],
	    dsubgrid_p_u[nDOF_trial_element],
	    dsubgrid_p_v[nDOF_trial_element],
	    dsubgrid_p_w[nDOF_trial_element],
	    dsubgrid_u_p[nDOF_trial_element],
	    dsubgrid_u_u[nDOF_trial_element],
	    dsubgrid_v_p[nDOF_trial_element],
	    dsubgrid_v_v[nDOF_trial_element],
	    dsubgrid_w_p[nDOF_trial_element],
	    dsubgrid_w_w[nDOF_trial_element],
	    tau_0=0.0,
	    tau_1=0.0;
          //
          //calculate solution and gradients at quadrature points
          //
	  p=0.0;u=0.0;v=0.0;w=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_p[I]=0.0;grad_u[I]=0.0;grad_v[I]=0.0;grad_w[I]=0.0;
	    }
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_j=eN*nDOF_trial_element+j;
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
              p += valFromDOF_c(p_dof[p_l2g[eN_j]],p_trial[eN_k_j]);
              u += valFromDOF_c(u_dof[vel_l2g[eN_j]],vel_trial[eN_k_j]);
              v += valFromDOF_c(v_dof[vel_l2g[eN_j]],vel_trial[eN_k_j]);
              w += valFromDOF_c(w_dof[vel_l2g[eN_j]],vel_trial[eN_k_j]);

	      for (int I=0;I<nSpace;I++)
		{
		  grad_p[I] += gradFromDOF_c(p_dof[p_l2g[eN_j]],p_grad_trial[eN_k_j_nSpace+I]);
		  grad_u[I] += gradFromDOF_c(u_dof[vel_l2g[eN_j]],vel_grad_trial[eN_k_j_nSpace+I]);
		  grad_v[I] += gradFromDOF_c(v_dof[vel_l2g[eN_j]],vel_grad_trial[eN_k_j_nSpace+I]);
		  grad_w[I] += gradFromDOF_c(w_dof[vel_l2g[eN_j]],vel_grad_trial[eN_k_j_nSpace+I]);
		}
	    }
          //
          //calculate pde coefficients at quadrature points
          //
 
          double rho,nu,mu,H_rho,d_rho,H_mu,d_mu;
          H_rho = smoothedHeaviside(eps_rho,phi[eN_k]);
          d_rho = smoothedDirac(eps_rho,phi[eN_k]);
          H_mu = smoothedHeaviside(eps_mu,phi[eN_k]);
          d_mu = smoothedDirac(eps_mu,phi[eN_k]);
  
          rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
          nu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
          mu  = rho_0*nu_0*(1.0-H_mu)+rho_1*nu_1*H_mu;
 
           //u momentum accumulation
          mom_u_acc=u;
          dmom_u_acc_u=1.0;
  
          //v momentum accumulation
          mom_v_acc=v;
          dmom_v_acc_v=1.0;
  
          //w momentum accumulation
          mom_w_acc=w;
          dmom_w_acc_w=1.0;          
          
          //
          //calculate time derivative at quadrature points
          //
          bdf_c(alphaBDF,
		q_mom_u_acc_beta_bdf[eN_k],
		mom_u_acc,
		dmom_u_acc_u,
		mom_u_acc_t,
		dmom_u_acc_u_t);
          bdf_c(alphaBDF,
		q_mom_v_acc_beta_bdf[eN_k],
		mom_v_acc,
		dmom_v_acc_v,
		mom_v_acc_t,
		dmom_v_acc_v_t);
          bdf_c(alphaBDF,
		q_mom_w_acc_beta_bdf[eN_k],
		mom_w_acc,
		dmom_w_acc_w,
		mom_w_acc_t,
		dmom_w_acc_w_t);          



         
          //
          //calculate subgrid error - strong residual
          //
          //calculate strong residual
	  pdeResidual_p = grad_u[0] + grad_v[1] + grad_w[2];

	  pdeResidual_u = rho *(mom_u_acc_t + u*grad_u[0] + v*grad_u[1] + w*grad_u[2] - g[0])
	                + grad_p[0]; // - mu poisson_u  ==> ommited as only first order right now !! 

	  pdeResidual_v = rho *(mom_v_acc_t + u*grad_v[0] + v*grad_v[1] + w*grad_v[2] - g[1])
	                + grad_p[1]; // - mu poisson_u  ==> ommited as only first order right now !! 
	                
	  pdeResidual_w = rho *(mom_w_acc_t + u*grad_w[0] + v*grad_w[1] + w*grad_w[2] - g[2])
	                + grad_p[2]; // - mu poisson_u  ==> ommited as only first order right now !! 
	                        

          //calculate tau and subgrid 
          calculateSubgridError_tau_c(hFactor,elementDiameter[eN],
				      dmom_u_acc_u_t,rho,mu,
				      u,v,w,
				      tau_0,tau_1,q_cfl[eN_k]);				      
				      

          subgrid_u = -tau_0*pdeResidual_u;
          subgrid_v = -tau_0*pdeResidual_v;
          subgrid_w = -tau_0*pdeResidual_w;

          subgrid_p = -tau_1*pdeResidual_p;

          
          //calculate the Jacobian of strong residual
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;

	      dpdeResidual_p_u[j]=vel_grad_trial[eN_k_j_nSpace+0];
	      dpdeResidual_p_v[j]=vel_grad_trial[eN_k_j_nSpace+1];
	      dpdeResidual_p_w[j]=vel_grad_trial[eN_k_j_nSpace+2];

	      dpdeResidual_u_p[j]=p_grad_trial[eN_k_j_nSpace+0];
	      dpdeResidual_u_u[j]=rho*(dmom_u_acc_u_t*vel_trial[eN_k_j]
	                             + u*vel_grad_trial[eN_k_j_nSpace+0] 
	                             + v*vel_grad_trial[eN_k_j_nSpace+1] 
	                             + w*vel_grad_trial[eN_k_j_nSpace+2]);

	      dpdeResidual_v_p[j]=p_grad_trial[eN_k_j_nSpace+1];
	      dpdeResidual_v_v[j]=rho*(dmom_v_acc_v_t*vel_trial[eN_k_j]
	                             + u*vel_grad_trial[eN_k_j_nSpace+0] 
	                             + v*vel_grad_trial[eN_k_j_nSpace+1] 
	                             + w*vel_grad_trial[eN_k_j_nSpace+2]);
	                             
	      dpdeResidual_w_p[j]=p_grad_trial[eN_k_j_nSpace+2];
	      dpdeResidual_w_w[j]=rho*(dmom_w_acc_w_t*vel_trial[eN_k_j]
	                             + u*vel_grad_trial[eN_k_j_nSpace+0] 
	                             + v*vel_grad_trial[eN_k_j_nSpace+1] 
	                             + w*vel_grad_trial[eN_k_j_nSpace+2]);
        /* p */
        dsubgrid_p_u[j] = -tau_1*dpdeResidual_p_u[j];
        dsubgrid_p_v[j] = -tau_1*dpdeResidual_p_v[j];
        dsubgrid_p_w[j] = -tau_1*dpdeResidual_p_w[j];

        /* u */
        dsubgrid_u_p[j] = -tau_0*dpdeResidual_u_p[j];
        dsubgrid_u_u[j] = -tau_0*dpdeResidual_u_u[j];
        /* v */
        dsubgrid_v_p[j] = -tau_0*dpdeResidual_v_p[j];
        dsubgrid_v_v[j] = -tau_0*dpdeResidual_v_v[j];
        /* w */
        dsubgrid_w_p[j] = -tau_0*dpdeResidual_w_p[j];
        dsubgrid_w_w[j] = -tau_0*dpdeResidual_w_w[j];
      }
    
    
    
    
  	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_k_i=eN_k*nDOF_test_element+i;
	      int eN_k_i_nSpace=eN_k_i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int eN_k_j=eN_k*nDOF_trial_element+j;
		  int eN_k_j_nSpace = eN_k_j*nSpace;

          /* p */
		  elementJacobian_p_p[i][j] -=  dsubgrid_u_p[j]*p_grad_test_dV[eN_k_i_nSpace+0]
	                                 +  dsubgrid_v_p[j]*p_grad_test_dV[eN_k_i_nSpace+1]
	                                 +  dsubgrid_w_p[j]*p_grad_test_dV[eN_k_i_nSpace+2];

		  elementJacobian_p_u[i][j] += vel_grad_trial[eN_k_j_nSpace+0]*p_test_dV[eN_k_i];
		  elementJacobian_p_v[i][j] += vel_grad_trial[eN_k_j_nSpace+1]*p_test_dV[eN_k_i];
		  elementJacobian_p_w[i][j] += vel_grad_trial[eN_k_j_nSpace+2]*p_test_dV[eN_k_i];

          /* u */
		  elementJacobian_u_p[i][j] -= p_trial[eN_k_j]*vel_grad_test_dV[eN_k_i_nSpace+0];

		  elementJacobian_u_u[i][j] += dpdeResidual_u_u[j]*vel_test_dV[eN_k_i]
	                            + mu*(vel_grad_trial[eN_k_j_nSpace+0]*vel_grad_test_dV[eN_k_i_nSpace+0]
	                                 +vel_grad_trial[eN_k_j_nSpace+1]*vel_grad_test_dV[eN_k_i_nSpace+1]
	                                 +vel_grad_trial[eN_k_j_nSpace+2]*vel_grad_test_dV[eN_k_i_nSpace+2])		 
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+0] * (dsubgrid_u_u[j]*u)// + vel_trial[eN_k_j]*subgrid_u + dsubgrid_u_u[j]*subgrid_u)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+1] * (dsubgrid_u_u[j]*v)// + vel_trial[eN_k_j]*subgrid_v + dsubgrid_u_u[j]*subgrid_v)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+2] * (dsubgrid_u_u[j]*w);// + vel_trial[eN_k_j]*subgrid_w + dsubgrid_u_u[j]*subgrid_w);
		  
		  elementJacobian_u_v[i][j] += 0.0;
		  elementJacobian_u_w[i][j] += 0.0;

          /* v */
		  elementJacobian_v_p[i][j] -= p_trial[eN_k_j]*vel_grad_test_dV[eN_k_i_nSpace+1];

		  elementJacobian_v_u[i][j] += 0.0; 		  
		  elementJacobian_v_v[i][j] += dpdeResidual_v_v[j]*vel_test_dV[eN_k_i]
	                            + mu*(vel_grad_trial[eN_k_j_nSpace+0]*vel_grad_test_dV[eN_k_i_nSpace+0]
	                                 +vel_grad_trial[eN_k_j_nSpace+1]*vel_grad_test_dV[eN_k_i_nSpace+1]
	                                 +vel_grad_trial[eN_k_j_nSpace+2]*vel_grad_test_dV[eN_k_i_nSpace+2])		 
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+0] * (dsubgrid_v_v[j]*u)// + vel_trial[eN_k_j]*subgrid_u + dsubgrid_v_v[j]*subgrid_u)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+1] * (dsubgrid_v_v[j]*v)// + vel_trial[eN_k_j]*subgrid_v + dsubgrid_v_v[j]*subgrid_v)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+2] * (dsubgrid_v_v[j]*w);// + vel_trial[eN_k_j]*subgrid_w + dsubgrid_v_v[j]*subgrid_w);
		  elementJacobian_v_w[i][j] += 0.0;
		 
          /* w */
		  elementJacobian_w_p[i][j] -= p_trial[eN_k_j]*vel_grad_test_dV[eN_k_i_nSpace+2];


		  elementJacobian_w_u[i][j] += 0.0;
		  elementJacobian_w_v[i][j] += 0.0;

		  elementJacobian_w_w[i][j] += dpdeResidual_w_w[j]*vel_test_dV[eN_k_i]
	                            + mu*(vel_grad_trial[eN_k_j_nSpace+0]*vel_grad_test_dV[eN_k_i_nSpace+0]
	                                 +vel_grad_trial[eN_k_j_nSpace+1]*vel_grad_test_dV[eN_k_i_nSpace+1]
	                                 +vel_grad_trial[eN_k_j_nSpace+2]*vel_grad_test_dV[eN_k_i_nSpace+2])		 
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+0] * (dsubgrid_w_w[j]*u)// + vel_trial[eN_k_j]*subgrid_u + dsubgrid_w_w[j]*subgrid_u)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+1] * (dsubgrid_w_w[j]*v)// + vel_trial[eN_k_j]*subgrid_v + dsubgrid_w_w[j]*subgrid_v)
		                        -rho*vel_grad_test_dV[eN_k_i_nSpace+2] * (dsubgrid_w_w[j]*w);// + vel_trial[eN_k_j]*subgrid_w + dsubgrid_w_w[j]*subgrid_w);
		  



		}//j
            }//i
	}//k
      //
      //load into element Jacobian into global Jacobian
      //
      for (int i=0;i<nDOF_test_element;i++)
	{
	  int eN_i = eN*nDOF_test_element+i;
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      int eN_i_j = eN_i*nDOF_trial_element+j;
	      globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_p_p[eN_i_j]] += elementJacobian_p_p[i][j];
	      globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_p_u[eN_i_j]] += elementJacobian_p_u[i][j];
	      globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_p_v[eN_i_j]] += elementJacobian_p_v[i][j];
	      globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_p_w[eN_i_j]] += elementJacobian_p_w[i][j];

	      globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_u_p[eN_i_j]] += elementJacobian_u_p[i][j];
	      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
	      globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_u_v[eN_i_j]] += elementJacobian_u_v[i][j];
	      globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_u_w[eN_i_j]] += elementJacobian_u_w[i][j];

	      globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_v_p[eN_i_j]] += elementJacobian_v_p[i][j];
	      globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_v_u[eN_i_j]] += elementJacobian_v_u[i][j];
	      globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_v_v[eN_i_j]] += elementJacobian_v_v[i][j];
	      globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_v_w[eN_i_j]] += elementJacobian_v_w[i][j];

	      globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_w_p[eN_i_j]] += elementJacobian_w_p[i][j];
	      globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_w_u[eN_i_j]] += elementJacobian_w_u[i][j];
	      globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_w_v[eN_i_j]] += elementJacobian_w_v[i][j];
	      globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_w_w[eN_i_j]] += elementJacobian_w_w[i][j];
	    }//j
	}//i
    }//elements
  //
  //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
  //
   for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      register int ebN = exteriorElementBoundariesArray[ebNE]; 
      register int eN  = elementBoundaryElementsArray[ebN*2+0];
      const double eps_rho = epsFact_rho*elementDiameter[eN],
      	eps_mu = epsFact_mu*elementDiameter[eN];
        register double  boundaryJacobian_p_p[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_p_u[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_p_v[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_p_w[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_u_p[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_u_u[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_u_v[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_u_w[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_v_p[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_v_u[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_v_v[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_v_w[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_w_p[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_w_u[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_w_v[nDOF_test_element][nDOF_trial_element],
	boundaryJacobian_w_w[nDOF_test_element][nDOF_trial_element];
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    boundaryJacobian_p_p[i][j]=0.0;
	    boundaryJacobian_p_u[i][j]=0.0;
	    boundaryJacobian_p_v[i][j]=0.0;
	    boundaryJacobian_p_w[i][j]=0.0;
	    boundaryJacobian_u_p[i][j]=0.0;
	    boundaryJacobian_u_u[i][j]=0.0;
	    boundaryJacobian_u_v[i][j]=0.0;
	    boundaryJacobian_u_w[i][j]=0.0;
	    boundaryJacobian_v_p[i][j]=0.0;
	    boundaryJacobian_v_u[i][j]=0.0;
	    boundaryJacobian_v_v[i][j]=0.0;
	    boundaryJacobian_v_w[i][j]=0.0;
	    boundaryJacobian_w_p[i][j]=0.0;
	    boundaryJacobian_w_u[i][j]=0.0;
	    boundaryJacobian_w_v[i][j]=0.0;
	    boundaryJacobian_w_w[i][j]=0.0;
	  }
  
  
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;

	  register double p_ext=0.0,
	    u_ext=0.0,
	    v_ext=0.0,
	    w_ext=0.0,
	    grad_p_ext[nSpace],
	    grad_u_ext[nSpace],
	    grad_v_ext[nSpace],
	    grad_w_ext[nSpace];
           double norm[nSpace];
	   	
	  norm[nSpace] = ebqe_n_ext[ebNE_kb_nSpace]; 	    
	    
	    
	  // 
	  //calculate the solution and gradients at quadrature points 
	  p_ext=0.0;u_ext=0.0;v_ext=0.0;w_ext=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_p_ext[I] = 0.0;
	      grad_u_ext[I] = 0.0;
	      grad_v_ext[I] = 0.0;
	      grad_w_ext[I] = 0.0;
	    }
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      register int eN_j = eN*nDOF_trial_element+j,
		ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
		ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	      p_ext += valFromDOF_c(p_dof[p_l2g[eN_j]],p_trial_ext[ebNE_kb_j]); 
	      u_ext += valFromDOF_c(u_dof[vel_l2g[eN_j]],vel_trial_ext[ebNE_kb_j]); 
	      v_ext += valFromDOF_c(v_dof[vel_l2g[eN_j]],vel_trial_ext[ebNE_kb_j]); 
	      w_ext += valFromDOF_c(w_dof[vel_l2g[eN_j]],vel_trial_ext[ebNE_kb_j]); 
               
	      for (int I=0;I<nSpace;I++)
		{
		  grad_p_ext[I] += gradFromDOF_c(p_dof[p_l2g[eN_j]],p_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_u_ext[I] += gradFromDOF_c(u_dof[vel_l2g[eN_j]],vel_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_v_ext[I] += gradFromDOF_c(v_dof[vel_l2g[eN_j]],vel_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_w_ext[I] += gradFromDOF_c(w_dof[vel_l2g[eN_j]],vel_grad_trial_ext[ebNE_kb_j_nSpace+I]);
		} 
	    }

          double rho,nu,mu,H_rho,d_rho,H_mu,d_mu,gamma;
          H_rho = smoothedHeaviside(eps_rho,phi[eN]);
          d_rho = smoothedDirac(eps_rho,phi[eN]);
          H_mu = smoothedHeaviside(eps_mu,phi[eN]);
          d_mu = smoothedDirac(eps_mu,phi[eN]);
  
          rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
          nu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
          mu  = rho_0*nu_0*(1.0-H_mu)+rho_1*nu_1*H_mu;
	            
          calculateInteriorPenalty   (hFactor,elementDiameter[eN],
				      mu,gamma);


	  //
	  //load the boundary values
	  //

         for (int i=0;i<nDOF_test_element;i++)
           {
             register int eN_i = eN*nDOF_test_element+i,
	       ebNE_kb_i = ebNE_kb*nDOF_test_element+i,
	       ebNE_kb_i_nSpace= ebNE_kb_i*nSpace; 
             for (int j=0;j<nDOF_trial_element;j++)
	       {
	         register int eN_j = eN*nDOF_trial_element+j,
		  ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
		  ebNE_kb_j_nSpace= ebNE_kb_j*nSpace; 
			       
	         boundaryJacobian_p_p[i][j]+= 0.0;
	         boundaryJacobian_p_u[i][j]+= vel_trial_ext[ebNE_kb_j]*norm[0]*p_test_dS_ext[ebNE_kb_i];
	         boundaryJacobian_p_v[i][j]+= vel_trial_ext[ebNE_kb_j]*norm[1]*p_test_dS_ext[ebNE_kb_i];
	         boundaryJacobian_p_w[i][j]+= vel_trial_ext[ebNE_kb_j]*norm[2]*p_test_dS_ext[ebNE_kb_i];
	     
	         boundaryJacobian_u_p[i][j]+= p_trial_ext[ebNE_kb_j]*norm[0]*vel_test_dS_ext[ebNE_kb_i];
	         boundaryJacobian_u_u[i][j]+=- mu*(vel_grad_trial_ext[ebNE_kb_j_nSpace+0]*norm[0] 
	                                         + vel_grad_trial_ext[ebNE_kb_j_nSpace+1]*norm[1] 
						 + vel_grad_trial_ext[ebNE_kb_j_nSpace+2]*norm[2])*vel_test_dS_ext[ebNE_kb_i]
					 //   - mu*vel_trial_ext[ebNE_kb_j]*(vel_grad_test_dS_ext[ebNE_kb_i_nSpace+0]*norm[0] 
  				         //                                + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+1]*norm[1] 
  				         //                                + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+2]*norm[2])
			                    + gamma*vel_trial_ext[ebNE_kb_j]*vel_test_dS_ext[ebNE_kb_i];		
	         boundaryJacobian_u_v[i][j]+= 0.0;
	         boundaryJacobian_u_w[i][j]+= 0.0;
		 
	     
	         boundaryJacobian_v_p[i][j]+= p_trial_ext[ebNE_kb_j]*norm[1]*vel_test_dS_ext[ebNE_kb_i];
	         boundaryJacobian_v_u[i][j]+= 0.0;
	         boundaryJacobian_v_v[i][j]+=- mu*(vel_grad_trial_ext[ebNE_kb_j_nSpace+0]*norm[0] 
	                                         + vel_grad_trial_ext[ebNE_kb_j_nSpace+1]*norm[1] 
						 + vel_grad_trial_ext[ebNE_kb_j_nSpace+2]*norm[2])*vel_test_dS_ext[ebNE_kb_i]
					//  - mu*vel_trial_ext[ebNE_kb_j]*(vel_grad_test_dS_ext[ebNE_kb_i_nSpace+0]*norm[0] 
  				        //				 + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+1]*norm[1] 
  				        //				 + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+2]*norm[2])
			                    + gamma*vel_trial_ext[ebNE_kb_j]*vel_test_dS_ext[ebNE_kb_i];
	         boundaryJacobian_v_w[i][j]+= 0.0;

	     
	         boundaryJacobian_w_p[i][j]+= p_trial_ext[ebNE_kb_j]*norm[2]*vel_test_dS_ext[ebNE_kb_i];
	         boundaryJacobian_w_u[i][j]+= 0.0;
		 boundaryJacobian_w_v[i][j]+= 0.0;
	         boundaryJacobian_w_w[i][j]+=- mu*(vel_grad_trial_ext[ebNE_kb_j_nSpace+0]*norm[0] 
	                                         + vel_grad_trial_ext[ebNE_kb_j_nSpace+1]*norm[1] 
						 + vel_grad_trial_ext[ebNE_kb_j_nSpace+2]*norm[2])*vel_test_dS_ext[ebNE_kb_i]
					 //  - mu*vel_trial_ext[ebNE_kb_j]*(vel_grad_test_dS_ext[ebNE_kb_i_nSpace+0]*norm[0] 
  				         //				  + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+1]*norm[1] 
  				         //				  + vel_grad_test_dS_ext[ebNE_kb_i_nSpace+2]*norm[2])
			                    + gamma*vel_trial_ext[ebNE_kb_j]*vel_test_dS_ext[ebNE_kb_i];
	         

	       }//j
           }//i	       // 


	  
	  
        }//kb
	  
       //
       //update the global Jacobian from the flux Jacobian
       //      
       for (int i=0;i<nDOF_test_element;i++)
         {
           register int eN_i = eN*nDOF_test_element+i;
           for (int j=0;j<nDOF_trial_element;j++)
	     {
	       register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
	       
	       globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_eb_p_p[ebN_i_j]] += boundaryJacobian_p_p[i][j];
	       globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_eb_p_u[ebN_i_j]] += boundaryJacobian_p_u[i][j];
	       globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_eb_p_v[ebN_i_j]] += boundaryJacobian_p_v[i][j];
	       globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_eb_p_w[ebN_i_j]] += boundaryJacobian_p_w[i][j];
	     	
	       globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_eb_u_p[ebN_i_j]] += boundaryJacobian_u_p[i][j];
	       globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += boundaryJacobian_u_u[i][j];
	       globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_eb_u_v[ebN_i_j]] += boundaryJacobian_u_v[i][j];
	       globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_eb_u_w[ebN_i_j]] += boundaryJacobian_u_w[i][j];
	     	
	       globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_eb_v_p[ebN_i_j]] += boundaryJacobian_v_p[i][j];
	       globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_eb_v_u[ebN_i_j]] += boundaryJacobian_v_u[i][j];
	       globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_eb_v_v[ebN_i_j]] += boundaryJacobian_v_v[i][j];
	       globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_eb_v_w[ebN_i_j]] += boundaryJacobian_v_w[i][j];
	     	
	       globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_eb_w_p[ebN_i_j]] += boundaryJacobian_w_p[i][j];
	       globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_eb_w_u[ebN_i_j]] += boundaryJacobian_w_u[i][j];
	       globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_eb_w_v[ebN_i_j]] += boundaryJacobian_w_v[i][j];
	       globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_eb_w_w[ebN_i_j]] += boundaryJacobian_w_w[i][j];

	     }//j
         }//i	       // 

	
    }//ebNE
    
}//computeJacobian

extern "C" void calculateVelocityAverage_RBLES2P(int nExteriorElementBoundaries_global,
						 int* exteriorElementBoundariesArray,
						 int nInteriorElementBoundaries_global,
						 int* interiorElementBoundariesArray,
						 int* elementBoundaryElementsArray,
						 int* elementBoundaryLocalElementBoundariesArray,
						 int* vel_l2g, 
						 double* u_dof, double* v_dof, double* w_dof,
						 double* vel_trial,
						 double* ebqe_velocity,
						 double* velocityAverage)
{
  using namespace RBLES2P;
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      register int ebN = exteriorElementBoundariesArray[ebNE],
	eN_global   = elementBoundaryElementsArray[ebN*2+0],
	ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace,
	    ebNE_kb_nSpace = ebNE*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
	  velocityAverage[ebN_kb_nSpace+0]=ebqe_velocity[ebNE_kb_nSpace+0];
	  velocityAverage[ebN_kb_nSpace+1]=ebqe_velocity[ebNE_kb_nSpace+1];
	  velocityAverage[ebN_kb_nSpace+2]=ebqe_velocity[ebNE_kb_nSpace+2];
	}//ebNE
    }
  for (int ebNI = 0; ebNI < nInteriorElementBoundaries_global; ebNI++) 
    { 
      register int ebN = interiorElementBoundariesArray[ebNI],
	left_eN_global   = elementBoundaryElementsArray[ebN*2+0],
	left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	right_eN_global  = elementBoundaryElementsArray[ebN*2+1],
	right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2+1];

      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
	  register double u_left=0.0,
	    v_left=0.0,
	    w_left=0.0,
	    u_right=0.0,
	    v_right=0.0,
	    w_right=0.0;
	  // 
	  //calculate the velocity solution at quadrature points on left and right
	  // 
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      int left_eN_j = left_eN_global*nDOF_trial_element+j;
	      int left_eN_ebN_kb_j = left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		kb*nDOF_trial_element +
		j;
	      int right_eN_j = right_eN_global*nDOF_trial_element+j;
	      int right_eN_ebN_kb_j = right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		kb*nDOF_trial_element +
		j;
	      u_left += valFromDOF_c(u_dof[vel_l2g[left_eN_j]],vel_trial[left_eN_ebN_kb_j]); 
	      v_left += valFromDOF_c(v_dof[vel_l2g[left_eN_j]],vel_trial[left_eN_ebN_kb_j]); 
	      w_left += valFromDOF_c(w_dof[vel_l2g[left_eN_j]],vel_trial[left_eN_ebN_kb_j]); 
	      u_right += valFromDOF_c(u_dof[vel_l2g[right_eN_j]],vel_trial[right_eN_ebN_kb_j]); 
	      v_right += valFromDOF_c(v_dof[vel_l2g[right_eN_j]],vel_trial[right_eN_ebN_kb_j]); 
	      w_right += valFromDOF_c(w_dof[vel_l2g[right_eN_j]],vel_trial[right_eN_ebN_kb_j]); 
	    }
	  velocityAverage[ebN_kb_nSpace+0]=0.5*(u_left + u_right);
	  velocityAverage[ebN_kb_nSpace+1]=0.5*(v_left + v_right);
	  velocityAverage[ebN_kb_nSpace+2]=0.5*(w_left + w_right);
	}//ebNI
    }






}
