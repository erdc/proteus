#include "RANS2PV2.h"
#include <iostream>
#include <cassert>

extern "C" void calculateResidual_RANS2PV2(//element
					   double* mesh_trial_ref,
					   double* mesh_grad_trial_ref,
					   double* mesh_dof,
					   int* mesh_l2g,
					   double* dV_ref,
					   double* p_trial_ref,
					   double* p_grad_trial_ref,
					   double* p_test_ref,
					   double* p_grad_test_ref,
					   double* vel_trial_ref,
					   double* vel_grad_trial_ref,
					   double* vel_test_ref,
					   double* vel_grad_test_ref,
					   //element boundary
					   double* mesh_trial_trace_ref,
					   double* mesh_grad_trial_trace_ref,
					   double* dS_ref,
					   double* p_trial_trace_ref,
					   double* p_grad_trial_trace_ref,
					   double* p_test_trace_ref,
					   double* p_grad_test_trace_ref,
					   double* vel_trial_trace_ref,
					   double* vel_grad_trial_trace_ref,
					   double* vel_test_trace_ref,
					   double* vel_grad_test_trace_ref,					 
					   double* normal_ref,
					   double* boundaryJac_ref,
					   //physics
					   int nElements_global,
					   double alphaBDF,
					   double epsFact_rho,
					   double epsFact_mu, 
					   double sigma,
					   double rho_0,
					   double nu_0,
					   double rho_1,
					   double nu_1,
					   double Ct_sge,
					   double Cd_sge,
					   double C_dc,
					   int* p_l2g, 
					   int* vel_l2g, 
					   double* p_dof, 
					   double* u_dof, 
					   double* v_dof, 
					   double* w_dof,
					   double* g,
					   double* phi,
					   double* normal_phi,
					   double* kappa_phi,
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
					   int offset_p, int offset_u, int offset_v, int offset_w, 
					   int stride_p, int stride_u, int stride_v, int stride_w, 
					   double* globalResidual,
					   int nExteriorElementBoundaries_global,
					   int* exteriorElementBoundariesArray,
					   int* elementBoundaryElementsArray,
					   int* elementBoundaryLocalElementBoundariesArray,
					   double* ebqe_phi_ext,
					   double* ebqe_normal_phi_ext,
					   double* ebqe_kappa_phi_ext,
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
					   double* q_velocity,
					   double* ebqe_velocity,
					   double* flux)
{
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;
  //
  //loop over elements to compute volume integrals and load them into element and global residual
  //
  double globalConservationError=0.0;
  for(int eN=0;eN<nElements_global;eN++)
    {
      //declare local storage for element residual and initialize
      register double elementResidual_p[nDOF_test_element],
	elementResidual_u[nDOF_test_element],
	elementResidual_v[nDOF_test_element],
	elementResidual_w[nDOF_test_element],
	eps_rho,eps_mu;
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_p[i]=0.0;
	  elementResidual_u[i]=0.0;
	  elementResidual_v[i]=0.0;
	  elementResidual_w[i]=0.0;
	}//i
      //
      //loop over quadrature points and compute integrands
      //
      for(int k=0;k<nQuadraturePoints_element;k++)
        {
	  //compute indices and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k,
	    eN_k_nSpace = eN_k*nSpace,
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double p=0.0,u=0.0,v=0.0,w=0.0,
	    grad_p[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
	    mom_u_acc=0.0,
	    dmom_u_acc_u=0.0,
	    mom_v_acc=0.0,
	    dmom_v_acc_v=0.0,
	    mom_w_acc=0.0,
	    dmom_w_acc_w=0.0,
	    mass_adv[nSpace],
	    dmass_adv_u[nSpace],
	    dmass_adv_v[nSpace],
	    dmass_adv_w[nSpace],
	    mom_u_adv[nSpace],
	    dmom_u_adv_u[nSpace],
	    dmom_u_adv_v[nSpace],
	    dmom_u_adv_w[nSpace],
	    mom_v_adv[nSpace],
	    dmom_v_adv_u[nSpace],
	    dmom_v_adv_v[nSpace],
	    dmom_v_adv_w[nSpace],
	    mom_w_adv[nSpace],
	    dmom_w_adv_u[nSpace],
	    dmom_w_adv_v[nSpace],
	    dmom_w_adv_w[nSpace],
	    mom_u_diff_ten[nSpace],
	    mom_v_diff_ten[nSpace],
	    mom_w_diff_ten[nSpace],
	    mom_uv_diff_ten[1],
	    mom_uw_diff_ten[1],
	    mom_vu_diff_ten[1],
	    mom_vw_diff_ten[1],
	    mom_wu_diff_ten[1],
	    mom_wv_diff_ten[1],
	    mom_u_source=0.0,
	    mom_v_source=0.0,
	    mom_w_source=0.0,
	    mom_u_ham=0.0,
	    dmom_u_ham_grad_p[nSpace],
	    mom_v_ham=0.0,
	    dmom_v_ham_grad_p[nSpace],
	    mom_w_ham=0.0,
	    dmom_w_ham_grad_p[nSpace],
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
	    Lstar_u_p[nDOF_test_element],
	    Lstar_v_p[nDOF_test_element],
	    Lstar_w_p[nDOF_test_element],
	    Lstar_u_u[nDOF_test_element],
	    Lstar_v_v[nDOF_test_element],
	    Lstar_w_w[nDOF_test_element],
	    Lstar_p_u[nDOF_test_element],
	    Lstar_p_v[nDOF_test_element],
	    Lstar_p_w[nDOF_test_element],
	    subgridError_p=0.0,
	    subgridError_u=0.0,
	    subgridError_v=0.0,
	    subgridError_w=0.0,
	    tau_p=0.0,
	    tau_v=0.0,
	    jac[nSpace*nSpace],
	    jacDet,
	    jacInv[nSpace*nSpace],
	    p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
	    p_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
	    p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
	    dV,x,y,z,
	    G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv,h_phi, vel[nSpace];
	  //get jacobian, etc for mapping reference element
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
	  ck.calculateG(jacInv,G,G_dd_G,tr_G);
	  ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);
	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;
	  //get the trial function gradients
	  ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
	  ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	  //get the solution
	  ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p);
	  ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],u);
	  ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],v);
	  ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w);
	  //get the solution gradients
	  ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p);
	  ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_u);
	  ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_v);
	  ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV;
	      vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
	      for (int I=0;I<nSpace;I++)
		{
		  p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		  vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		}
	    }
	  //save velocity at quadrature points for other models to use
	  q_velocity[eN_k_nSpace+0]=u;
	  q_velocity[eN_k_nSpace+1]=v;
	  q_velocity[eN_k_nSpace+2]=w;
	  //std::cout<<"velocity "<<u<<'\t'<<v<<'\t'<<w<<std::endl;
	  // q_velocity_last[eN_k_nSpace+0]=1.0;
	  // q_velocity_last[eN_k_nSpace+1]=1.0;
	  // q_velocity_last[eN_k_nSpace+2]=1.0;
	  // q_velocity_last[eN_k_nSpace+0]=u;
	  // q_velocity_last[eN_k_nSpace+1]=v;
	  // q_velocity_last[eN_k_nSpace+2]=w;
	  // //
	  // //debugging section for finite element calculations on interior
	  // //
	  // std::cout<<"eN = "<<eN<<" k = "<<k<<std::endl;
	  // for (int j=0;j<nDOF_trial_element;j++)
	  //   {
	  //     std::cout<<"p_trial["<<j<<"] "<<p_trial_ref[k*nDOF_trial_element+j]<<std::endl;
	  //     std::cout<<"vel_trial["<<j<<"] "<<vel_trial_ref[k*nDOF_trial_element+j]<<std::endl;
	  //     for (int I=0;I<nSpace;I++)
	  // 	{
	  // 	  std::cout<<"p_grad_trial["<<j<<","<<I<<"] "<<p_grad_trial[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"vel_grad_trial["<<j<<","<<I<<"] "<<vel_grad_trial[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"p_grad_test_dV["<<j<<","<<I<<"] "<<p_grad_test_dV[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"vel_grad_test_dV["<<j<<","<<I<<"] "<<vel_grad_test_dV[j*nSpace+I]<<std::endl;
	  // 	}
	  //   }
	  // std::cout<<"p "<<p<<std::endl;
	  // std::cout<<"u "<<u<<std::endl;
	  // std::cout<<"v "<<v<<std::endl;
	  // std::cout<<"w "<<w<<std::endl;
	  // for(int I=0;I<nSpace;I++)
	  //   {
	  //     std::cout<<"grad_p["<<I<<"] "<<grad_p[I]<<std::endl;
	  //     std::cout<<"grad_u["<<I<<"] "<<grad_u[I]<<std::endl;
	  //     std::cout<<"grad_v["<<I<<"] "<<grad_v[I]<<std::endl;
	  //     std::cout<<"grad_w["<<I<<"] "<<grad_w[I]<<std::endl;
	  //   }
          //
          //calculate pde coefficients at quadrature points
          //
	  //cek debug
	  //eps_rho=0.1;
	  //eps_mu=0.1;
	  RANS2PV2::evaluateCoefficients(eps_rho,
					 eps_mu,
					 sigma,
					 rho_0,
					 nu_0,
					 rho_1,
					 nu_1,
					 g,
					 phi[eN_k],
					 &normal_phi[eN_k_nSpace],
					 kappa_phi[eN_k],
					 p,
					 grad_p,
					 u,
					 v,
					 w,
					 mom_u_acc,
					 dmom_u_acc_u,
					 mom_v_acc,
					 dmom_v_acc_v,
					 mom_w_acc,
					 dmom_w_acc_w,
					 mass_adv,
					 dmass_adv_u,
					 dmass_adv_v,
					 dmass_adv_w,
					 mom_u_adv,
				 	 dmom_u_adv_u,
					 dmom_u_adv_v,
					 dmom_u_adv_w,
					 mom_v_adv,
					 dmom_v_adv_u,
					 dmom_v_adv_v,
					 dmom_v_adv_w,
					 mom_w_adv,
					 dmom_w_adv_u,
					 dmom_w_adv_v,
					 dmom_w_adv_w,
					 mom_u_diff_ten,
					 mom_v_diff_ten,
					 mom_w_diff_ten,
					 mom_uv_diff_ten,
					 mom_uw_diff_ten,
					 mom_vu_diff_ten,
					 mom_vw_diff_ten,
					 mom_wu_diff_ten,
					 mom_wv_diff_ten,
					 mom_u_source,
					 mom_v_source,
					 mom_w_source,
					 mom_u_ham,
					 dmom_u_ham_grad_p,
					 mom_v_ham,
					 dmom_v_ham_grad_p,
					 mom_w_ham,
					 dmom_w_ham_grad_p);          
	  //
	  //save momentum for time history and velocity for subgrid error
	  //
	  q_mom_u_acc[eN_k] = mom_u_acc;                            
	  q_mom_v_acc[eN_k] = mom_v_acc;                            
	  q_mom_w_acc[eN_k] = mom_w_acc;
	  //subgrid error uses grid scale velocity
	  q_mass_adv[eN_k_nSpace+0] = u;
	  q_mass_adv[eN_k_nSpace+1] = v;
	  q_mass_adv[eN_k_nSpace+2] = w;
          //
          //moving mesh
          //
          //omit for now
          //
          //calculate time derivative at quadrature points
          //
          ck.bdf(alphaBDF,
		 q_mom_u_acc_beta_bdf[eN_k],
		 mom_u_acc,
		 dmom_u_acc_u,
		 mom_u_acc_t,
		 dmom_u_acc_u_t);
          ck.bdf(alphaBDF,
		 q_mom_v_acc_beta_bdf[eN_k],
		 mom_v_acc,
		 dmom_v_acc_v,
		 mom_v_acc_t,
		 dmom_v_acc_v_t);
          ck.bdf(alphaBDF,
		 q_mom_w_acc_beta_bdf[eN_k],
		 mom_w_acc,
		 dmom_w_acc_w,
		 mom_w_acc_t,
		 dmom_w_acc_w_t);
          //
          //calculate subgrid error (strong residual and adjoint)
          //
          //calculate strong residual
	  pdeResidual_p = ck.Advection_strong(dmass_adv_u,grad_u) +
	    ck.Advection_strong(dmass_adv_v,grad_v) +
	    ck.Advection_strong(dmass_adv_w,grad_w);
	  
	  pdeResidual_u = ck.Mass_strong(mom_u_acc_t) +
	    ck.Advection_strong(&q_velocity_last[eN_k_nSpace],grad_u) +
	    ck.Hamiltonian_strong(dmom_u_ham_grad_p,grad_p) +
	    ck.Reaction_strong(mom_u_source);
	  
	  pdeResidual_v = ck.Mass_strong(mom_v_acc_t) +
	    ck.Advection_strong(&q_velocity_last[eN_k_nSpace],grad_v) +
	    ck.Hamiltonian_strong(dmom_v_ham_grad_p,grad_p) + 
	    ck.Reaction_strong(mom_v_source);
	  
	  pdeResidual_w = ck.Mass_strong(mom_w_acc_t) + 
	    ck.Advection_strong(&q_velocity_last[eN_k_nSpace],grad_w) +
	    ck.Hamiltonian_strong(dmom_w_ham_grad_p,grad_p) +
	    ck.Reaction_strong(mom_w_source);
	
          //calculate tau and tau*Res
	    RANS2PV2::calculateSubgridError_tau(Ct_sge,
					      Cd_sge,
					      G,
					      G_dd_G,
					      tr_G,
					      dmom_u_acc_u,//rho
					      dmom_u_acc_u_t/dmom_u_acc_u,//Dt
	  				      &q_velocity_last[eN_k_nSpace],//v
					      mom_u_diff_ten[1],//mu
	  				      tau_p,
					      tau_v,
					      q_cfl[eN_k]);

        RANS2PV2::calculateSubgridError_tauRes(tau_p,
	  					 tau_v,
	  					 pdeResidual_p,
	  					 pdeResidual_u,
	  					 pdeResidual_v,
	  					 pdeResidual_w,
	  					 subgridError_p,
	  					 subgridError_u,
	  					 subgridError_v,
	  					 subgridError_w);
     	//cek/ido todo use add options for lagging velocity and turning of RBLES terms
	  
         // Velocity  ==> SUPG   
	    //vel[0] = u;
	    //vel[1] = v;
	    //vel[2] = w;
 
        // Velocity + subgrid ==> RBLES
	    vel[0] = u + subgridError_u;
	    vel[1] = v + subgridError_v;
	    vel[2] = w + subgridError_w;

        // adjoint times the test functions 
          for (int i=0;i<nDOF_test_element;i++)
            {
	      register int i_nSpace = i*nSpace;
	      Lstar_u_p[i]=ck.Advection_adjoint(dmass_adv_u,&p_grad_test_dV[i_nSpace]);
	      Lstar_v_p[i]=ck.Advection_adjoint(dmass_adv_v,&p_grad_test_dV[i_nSpace]);
	      Lstar_w_p[i]=ck.Advection_adjoint(dmass_adv_w,&p_grad_test_dV[i_nSpace]);
	      Lstar_u_u[i]=ck.Advection_adjoint(vel,&vel_grad_test_dV[i_nSpace]);
	      Lstar_v_v[i]=ck.Advection_adjoint(vel,&vel_grad_test_dV[i_nSpace]);
	      Lstar_w_w[i]=ck.Advection_adjoint(vel,&vel_grad_test_dV[i_nSpace]);
	      Lstar_p_u[i]=ck.Hamiltonian_adjoint(dmom_u_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
	      Lstar_p_v[i]=ck.Hamiltonian_adjoint(dmom_v_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
	      Lstar_p_w[i]=ck.Hamiltonian_adjoint(dmom_w_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
            }

	  //cek todo
	  // mom_u_adv[0] = (u + subgridError_u)*(u + subgridError_u)
	  // mom_u_adv[1] = (u + subgridError_u)*(v + subgridError_v)
	  // mom_u_adv[2] = (u + subgridError_u)*(w + subgridError_w)
	  //and so on for mom_v_adv and mom_w_adv
	  //BUT THEN FIX SUBGRIDERROR or add a new term (e.g. RBLES_TERM)
	  //end todo
          //calculate shock capturing diffusion
          
	  norm_Rv = sqrt(pdeResidual_u*pdeResidual_u + pdeResidual_v*pdeResidual_v + pdeResidual_w*pdeResidual_w);
	  q_numDiff_u[eN_k] = C_dc*norm_Rv/sqrt(G_dd_G);
	  q_numDiff_v[eN_k] = q_numDiff_u[eN_k];
	  q_numDiff_w[eN_k] = q_numDiff_u[eN_k];
	  // //cek debug
	  // q_numDiff_u[eN_k] = 0.0;
	  // q_numDiff_v[eN_k] = 0.0;
	  // q_numDiff_w[eN_k] = 0.0;
	  // q_numDiff_u_last[eN_k] = 0.0;	  
	  // q_numDiff_v_last[eN_k] = 0.0;	  
	  // q_numDiff_w_last[eN_k] = 0.0;	  
          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int i_nSpace=i*nSpace;

	      elementResidual_p[i] += ck.Advection_weak(mass_adv,&p_grad_test_dV[i_nSpace]) +
		ck.SubgridError(subgridError_u,Lstar_u_p[i]) + 
		ck.SubgridError(subgridError_v,Lstar_v_p[i]) + 
		ck.SubgridError(subgridError_w,Lstar_w_p[i]);

	      elementResidual_u[i] += ck.Mass_weak(mom_u_acc_t,vel_test_dV[i]) + 
		ck.Advection_weak(mom_u_adv,&vel_grad_test_dV[i_nSpace]) +
		ck.Diffusion_weak(sdInfo_u_u_rowptr,sdInfo_u_u_colind,mom_u_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) + 
		ck.Diffusion_weak(sdInfo_u_v_rowptr,sdInfo_u_v_colind,mom_uv_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) + 
		ck.Diffusion_weak(sdInfo_u_w_rowptr,sdInfo_u_w_colind,mom_uw_diff_ten,grad_w,&vel_grad_test_dV[i_nSpace]) + 
		ck.Reaction_weak(mom_u_source,vel_test_dV[i]) + 
		ck.Hamiltonian_weak(mom_u_ham,vel_test_dV[i]) + 
		ck.SubgridError(subgridError_p,Lstar_p_u[i]) + 
		ck.SubgridError(subgridError_u,Lstar_u_u[i]) + 
		ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u,&vel_grad_test_dV[i_nSpace]); 
		 
	      elementResidual_v[i] += ck.Mass_weak(mom_v_acc_t,vel_test_dV[i]) + 
		ck.Advection_weak(mom_v_adv,&vel_grad_test_dV[i_nSpace]) +
		ck.Diffusion_weak(sdInfo_v_v_rowptr,sdInfo_v_v_colind,mom_v_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) + 
		ck.Diffusion_weak(sdInfo_v_u_rowptr,sdInfo_v_u_colind,mom_vu_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) + 
		ck.Diffusion_weak(sdInfo_v_w_rowptr,sdInfo_v_w_colind,mom_vw_diff_ten,grad_w,&vel_grad_test_dV[i_nSpace]) + 
		ck.Reaction_weak(mom_v_source,vel_test_dV[i]) + 
		ck.Hamiltonian_weak(mom_v_ham,vel_test_dV[i]) + 
		ck.SubgridError(subgridError_p,Lstar_p_v[i]) + 
		ck.SubgridError(subgridError_v,Lstar_v_v[i]) + 
		ck.NumericalDiffusion(q_numDiff_v_last[eN_k],grad_v,&vel_grad_test_dV[i_nSpace]); 

	      elementResidual_w[i] +=  ck.Mass_weak(mom_w_acc_t,vel_test_dV[i]) +
		ck.Advection_weak(mom_w_adv,&vel_grad_test_dV[i_nSpace]) + 
		ck.Diffusion_weak(sdInfo_w_w_rowptr,sdInfo_w_w_colind,mom_w_diff_ten,grad_w,&vel_grad_test_dV[i_nSpace]) + 
		ck.Diffusion_weak(sdInfo_w_u_rowptr,sdInfo_w_u_colind,mom_wu_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) + 
		ck.Diffusion_weak(sdInfo_w_v_rowptr,sdInfo_w_v_colind,mom_wv_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) + 
		ck.Reaction_weak(mom_w_source,vel_test_dV[i]) + 
		ck.Hamiltonian_weak(mom_w_ham,vel_test_dV[i]) + 
		ck.SubgridError(subgridError_p,Lstar_p_w[i]) + 
		ck.SubgridError(subgridError_w,Lstar_w_w[i]) + 
		ck.NumericalDiffusion(q_numDiff_w_last[eN_k],grad_w,&vel_grad_test_dV[i_nSpace]); 
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

          globalResidual[offset_p+stride_p*p_l2g[eN_i]]+=elementResidual_p[i];
          globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
          globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
          globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
        }//i
      // //
      // //debug
      // //
      // for(int i=0;i<nDOF_test_element;i++) 
      //   { 
      // 	  std::cout<<"eN "<<eN<<" i "<<i<<std::endl;
      // 	  std::cout<<"r_p"<<elementResidual_p[i]<<std::endl;
      // 	  std::cout<<"r_u"<<elementResidual_u[i]<<std::endl;
      // 	  std::cout<<"r_v"<<elementResidual_v[i]<<std::endl;
      // 	  std::cout<<"r_w"<<elementResidual_w[i]<<std::endl;
      // 	}
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
      register double elementResidual_p[nDOF_test_element],
	elementResidual_u[nDOF_test_element],
	elementResidual_v[nDOF_test_element],
	elementResidual_w[nDOF_test_element],
	eps_rho,eps_mu;
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
	    ebNE_kb_nSpace = ebNE_kb*nSpace,
	    ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
	    ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	  register double p_ext=0.0,
	    u_ext=0.0,
	    v_ext=0.0,
	    w_ext=0.0,
	    grad_p_ext[nSpace],
	    grad_u_ext[nSpace],
	    grad_v_ext[nSpace],
	    grad_w_ext[nSpace],
	    mom_u_acc_ext=0.0,
	    dmom_u_acc_u_ext=0.0,
	    mom_v_acc_ext=0.0,
	    dmom_v_acc_v_ext=0.0,
	    mom_w_acc_ext=0.0,
	    dmom_w_acc_w_ext=0.0,
	    mass_adv_ext[nSpace],
	    dmass_adv_u_ext[nSpace],
	    dmass_adv_v_ext[nSpace],
	    dmass_adv_w_ext[nSpace],
	    mom_u_adv_ext[nSpace],
	    dmom_u_adv_u_ext[nSpace],
	    dmom_u_adv_v_ext[nSpace],
	    dmom_u_adv_w_ext[nSpace],
	    mom_v_adv_ext[nSpace],
	    dmom_v_adv_u_ext[nSpace],
	    dmom_v_adv_v_ext[nSpace],
	    dmom_v_adv_w_ext[nSpace],
	    mom_w_adv_ext[nSpace],
	    dmom_w_adv_u_ext[nSpace],
	    dmom_w_adv_v_ext[nSpace],
	    dmom_w_adv_w_ext[nSpace],
	    mom_u_diff_ten_ext[nSpace],
	    mom_v_diff_ten_ext[nSpace],
	    mom_w_diff_ten_ext[nSpace],
	    mom_uv_diff_ten_ext[1],
	    mom_uw_diff_ten_ext[1],
	    mom_vu_diff_ten_ext[1],
	    mom_vw_diff_ten_ext[1],
	    mom_wu_diff_ten_ext[1],
	    mom_wv_diff_ten_ext[1],
	    mom_u_source_ext=0.0,
	    mom_v_source_ext=0.0,
	    mom_w_source_ext=0.0,
	    mom_u_ham_ext=0.0,
	    dmom_u_ham_grad_p_ext[nSpace],
	    mom_v_ham_ext=0.0,
	    dmom_v_ham_grad_p_ext[nSpace],
	    mom_w_ham_ext=0.0,
	    dmom_w_ham_grad_p_ext[nSpace],
	    dmom_u_adv_p_ext[nSpace],
	    dmom_v_adv_p_ext[nSpace],
	    dmom_w_adv_p_ext[nSpace],
	    flux_mass_ext=0.0,
	    flux_mom_u_adv_ext=0.0,
	    flux_mom_v_adv_ext=0.0,
	    flux_mom_w_adv_ext=0.0,
	    flux_mom_u_diff_ext=0.0,
	    flux_mom_v_diff_ext=0.0,
	    flux_mom_w_diff_ext=0.0,
	    bc_p_ext=0.0,
	    bc_u_ext=0.0,
	    bc_v_ext=0.0,
	    bc_w_ext=0.0,
	    bc_mom_u_acc_ext=0.0,
	    bc_dmom_u_acc_u_ext=0.0,
	    bc_mom_v_acc_ext=0.0,
	    bc_dmom_v_acc_v_ext=0.0,
	    bc_mom_w_acc_ext=0.0,
	    bc_dmom_w_acc_w_ext=0.0,
	    bc_mass_adv_ext[nSpace],
	    bc_dmass_adv_u_ext[nSpace],
	    bc_dmass_adv_v_ext[nSpace],
	    bc_dmass_adv_w_ext[nSpace],
	    bc_mom_u_adv_ext[nSpace],
	    bc_dmom_u_adv_u_ext[nSpace],
	    bc_dmom_u_adv_v_ext[nSpace],
	    bc_dmom_u_adv_w_ext[nSpace],
	    bc_mom_v_adv_ext[nSpace],
	    bc_dmom_v_adv_u_ext[nSpace],
	    bc_dmom_v_adv_v_ext[nSpace],
	    bc_dmom_v_adv_w_ext[nSpace],
	    bc_mom_w_adv_ext[nSpace],
	    bc_dmom_w_adv_u_ext[nSpace],
	    bc_dmom_w_adv_v_ext[nSpace],
	    bc_dmom_w_adv_w_ext[nSpace],
	    bc_mom_u_diff_ten_ext[nSpace],
	    bc_mom_v_diff_ten_ext[nSpace],
	    bc_mom_w_diff_ten_ext[nSpace],
	    bc_mom_uv_diff_ten_ext[1],
	    bc_mom_uw_diff_ten_ext[1],
	    bc_mom_vu_diff_ten_ext[1],
	    bc_mom_vw_diff_ten_ext[1],
	    bc_mom_wu_diff_ten_ext[1],
	    bc_mom_wv_diff_ten_ext[1],
	    bc_mom_u_source_ext=0.0,
	    bc_mom_v_source_ext=0.0,
	    bc_mom_w_source_ext=0.0,
	    bc_mom_u_ham_ext=0.0,
	    bc_dmom_u_ham_grad_p_ext[nSpace],
	    bc_mom_v_ham_ext=0.0,
	    bc_dmom_v_ham_grad_p_ext[nSpace],
	    bc_mom_w_ham_ext=0.0,
	    bc_dmom_w_ham_grad_p_ext[nSpace],
	    jac_ext[nSpace*nSpace],
	    jacDet_ext,
	    jacInv_ext[nSpace*nSpace],
	    boundaryJac[nSpace*(nSpace-1)],
	    metricTensor[(nSpace-1)*(nSpace-1)],
	    metricTensorDetSqrt,
	    dS,p_test_dS[nDOF_test_element],vel_test_dS[nDOF_test_element],
	    p_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_trial_element*nSpace],
	    normal[3],x_ext,y_ext,z_ext,
	    G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty;
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
	  dS = metricTensorDetSqrt*dS_ref[kb];
	  //get the metric tensor
	  //cek todo use symmetry
	  ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	  ck.calculateGScale(G,&ebqe_normal_phi_ext[ebNE_kb_nSpace],h_phi);
	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;
	  //compute shape and solution information
	  //shape
	  ck.gradTrialFromRef(&p_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
	  ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace);
	  //solution and gradients	
	  ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],p_ext);
	  ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	  ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext);
	  ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],w_ext);
	  ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
	  ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_u_ext);
	  ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_v_ext);
	  ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_w_ext);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      p_test_dS[j] = p_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	      vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	    }
	  // //
	  // //debugging section for finite element calculations on exterior
	  // //
	  // std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl;
	  // for (int j=0;j<nDOF_trial_element;j++)
	  //   {
	  //     std::cout<<"p_trial_trace["<<j<<"] "<<p_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	  //     std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	  //     std::cout<<"p_test_dS["<<j<<"] "<<p_test_dS[j]<<std::endl;
	  //     std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl;
	  //     for (int I=0;I<nSpace;I++)
	  // 	{
	  // 	  std::cout<<"p_grad_trial_trace["<<j<<","<<I<<"] "<<p_grad_trial_trace[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl;
	  // 	}
	  //   }
	  // std::cout<<"p_ext "<<p_ext<<std::endl;
	  // std::cout<<"u_ext "<<u_ext<<std::endl;
	  // std::cout<<"v_ext "<<v_ext<<std::endl;
	  // std::cout<<"w_ext "<<w_ext<<std::endl;
	  // for(int I=0;I<nSpace;I++)
	  //   {
	  //     std::cout<<"grad_p_ext["<<I<<"] "<<grad_p_ext[I]<<std::endl;
	  //     std::cout<<"grad_u_ext["<<I<<"] "<<grad_u_ext[I]<<std::endl;
	  //     std::cout<<"grad_v_ext["<<I<<"] "<<grad_v_ext[I]<<std::endl;
	  //     std::cout<<"grad_w_ext["<<I<<"] "<<grad_w_ext[I]<<std::endl;
	  //   }
	  //
	  //load the boundary values
	  //
	  bc_p_ext = isDOFBoundary_p[ebNE_kb]*ebqe_bc_p_ext[ebNE_kb]+(1-isDOFBoundary_p[ebNE_kb])*p_ext;
	  bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	  bc_v_ext = isDOFBoundary_v[ebNE_kb]*ebqe_bc_v_ext[ebNE_kb]+(1-isDOFBoundary_v[ebNE_kb])*v_ext;
	  bc_w_ext = isDOFBoundary_w[ebNE_kb]*ebqe_bc_w_ext[ebNE_kb]+(1-isDOFBoundary_w[ebNE_kb])*w_ext;
	  // 
	  //calculate the pde coefficients using the solution and the boundary values for the solution 
	  // 
	  //cek debug
	  //eps_rho=0.1;
	  //eps_mu=0.1;
	  RANS2PV2::evaluateCoefficients(eps_rho,
					 eps_mu,
					 sigma,
					 rho_0,
					 nu_0,
					 rho_1,
					 nu_1,
					 g,
					 ebqe_phi_ext[ebNE_kb],
					 &ebqe_normal_phi_ext[ebNE_kb_nSpace],
					 ebqe_kappa_phi_ext[ebNE_kb],
					 p_ext,
					 grad_p_ext,
					 u_ext,
					 v_ext,
					 w_ext,
					 mom_u_acc_ext,
					 dmom_u_acc_u_ext,
					 mom_v_acc_ext,
					 dmom_v_acc_v_ext,
					 mom_w_acc_ext,
					 dmom_w_acc_w_ext,
					 mass_adv_ext,
					 dmass_adv_u_ext,
					 dmass_adv_v_ext,
					 dmass_adv_w_ext,
					 mom_u_adv_ext,
					 dmom_u_adv_u_ext,
					 dmom_u_adv_v_ext,
					 dmom_u_adv_w_ext,
					 mom_v_adv_ext,
					 dmom_v_adv_u_ext,
					 dmom_v_adv_v_ext,
					 dmom_v_adv_w_ext,
					 mom_w_adv_ext,
					 dmom_w_adv_u_ext,
					 dmom_w_adv_v_ext,
					 dmom_w_adv_w_ext,
					 mom_u_diff_ten_ext,
					 mom_v_diff_ten_ext,
					 mom_w_diff_ten_ext,
					 mom_uv_diff_ten_ext,
					 mom_uw_diff_ten_ext,
					 mom_vu_diff_ten_ext,
					 mom_vw_diff_ten_ext,
					 mom_wu_diff_ten_ext,
					 mom_wv_diff_ten_ext,
					 mom_u_source_ext,
					 mom_v_source_ext,
					 mom_w_source_ext,
					 mom_u_ham_ext,
					 dmom_u_ham_grad_p_ext,
					 mom_v_ham_ext,
					 dmom_v_ham_grad_p_ext,
					 mom_w_ham_ext,
					 dmom_w_ham_grad_p_ext);          
	  RANS2PV2::evaluateCoefficients(eps_rho,
					 eps_mu,
					 sigma,
					 rho_0,
					 nu_0,
					 rho_1,
					 nu_1,
					 g,
					 ebqe_phi_ext[ebNE_kb],
					 &ebqe_normal_phi_ext[ebNE_kb_nSpace],
					 ebqe_kappa_phi_ext[ebNE_kb],
					 bc_p_ext,
					 grad_p_ext,//cek should't be used
					 bc_u_ext,
					 bc_v_ext,
					 bc_w_ext,
					 bc_mom_u_acc_ext,
					 bc_dmom_u_acc_u_ext,
					 bc_mom_v_acc_ext,
					 bc_dmom_v_acc_v_ext,
					 bc_mom_w_acc_ext,
					 bc_dmom_w_acc_w_ext,
					 bc_mass_adv_ext,
					 bc_dmass_adv_u_ext,
					 bc_dmass_adv_v_ext,
					 bc_dmass_adv_w_ext,
					 bc_mom_u_adv_ext,
					 bc_dmom_u_adv_u_ext,
					 bc_dmom_u_adv_v_ext,
					 bc_dmom_u_adv_w_ext,
					 bc_mom_v_adv_ext,
					 bc_dmom_v_adv_u_ext,
					 bc_dmom_v_adv_v_ext,
					 bc_dmom_v_adv_w_ext,
					 bc_mom_w_adv_ext,
					 bc_dmom_w_adv_u_ext,
					 bc_dmom_w_adv_v_ext,
					 bc_dmom_w_adv_w_ext,
					 bc_mom_u_diff_ten_ext,
					 bc_mom_v_diff_ten_ext,
					 bc_mom_w_diff_ten_ext,
					 bc_mom_uv_diff_ten_ext,
					 bc_mom_uw_diff_ten_ext,
					 bc_mom_vu_diff_ten_ext,
					 bc_mom_vw_diff_ten_ext,
					 bc_mom_wu_diff_ten_ext,
					 bc_mom_wv_diff_ten_ext,
					 bc_mom_u_source_ext,
					 bc_mom_v_source_ext,
					 bc_mom_w_source_ext,
					 bc_mom_u_ham_ext,
					 bc_dmom_u_ham_grad_p_ext,
					 bc_mom_v_ham_ext,
					 bc_dmom_v_ham_grad_p_ext,
					 bc_mom_w_ham_ext,
					 bc_dmom_w_ham_grad_p_ext);          
	  // 
	  //calculate the numerical fluxes 
	  // 
	  //cek debug
	  //ebqe_penalty_ext[ebNE_kb] = 10.0;
	  //
	  ck.calculateGScale(G,normal,h_penalty);
	  h_penalty = 10.0/h_penalty;
	  //cek debug
	  //h_penalty = 10.0/0.1;//he=0.1
	  RANS2PV2::exteriorNumericalAdvectiveFlux(isDOFBoundary_p[ebNE_kb],
						   isDOFBoundary_u[ebNE_kb],
						   isDOFBoundary_v[ebNE_kb],
						   isDOFBoundary_w[ebNE_kb],
						   isAdvectiveFluxBoundary_p[ebNE_kb],
						   isAdvectiveFluxBoundary_u[ebNE_kb],
						   isAdvectiveFluxBoundary_v[ebNE_kb],
						   isAdvectiveFluxBoundary_w[ebNE_kb],
						   normal,
						   bc_p_ext,
						   bc_mass_adv_ext,
						   bc_mom_u_adv_ext,
						   bc_mom_v_adv_ext,
						   bc_mom_w_adv_ext,
						   ebqe_bc_flux_mass_ext[ebNE_kb],
						   ebqe_bc_flux_mom_u_adv_ext[ebNE_kb],
						   ebqe_bc_flux_mom_v_adv_ext[ebNE_kb],
						   ebqe_bc_flux_mom_w_adv_ext[ebNE_kb],
						   p_ext,
						   mass_adv_ext,
						   mom_u_adv_ext,
						   mom_v_adv_ext,
						   mom_w_adv_ext,
						   dmass_adv_u_ext,
						   dmass_adv_v_ext,
						   dmass_adv_w_ext,
						   dmom_u_adv_p_ext,
						   dmom_u_adv_u_ext,
						   dmom_u_adv_v_ext,
						   dmom_u_adv_w_ext,
						   dmom_v_adv_p_ext,
						   dmom_v_adv_u_ext,
						   dmom_v_adv_v_ext,
						   dmom_v_adv_w_ext,
						   dmom_w_adv_p_ext,
						   dmom_w_adv_u_ext,
						   dmom_w_adv_v_ext,
						   dmom_w_adv_w_ext,
						   flux_mass_ext,
						   flux_mom_u_adv_ext,
						   flux_mom_v_adv_ext,
						   flux_mom_w_adv_ext,
						   &ebqe_velocity[ebNE_kb_nSpace]);
	  //cek todo need to switch to full stress and add adjoint consistency
	  RANS2PV2::exteriorNumericalDiffusiveFlux(eps_rho,
						   ebqe_phi_ext[ebNE_kb],
						   sdInfo_u_u_rowptr,
						   sdInfo_u_u_colind,
						   isDOFBoundary_u[ebNE_kb],
						   isDiffusiveFluxBoundary_u[ebNE_kb],
						   normal,
						   bc_mom_u_diff_ten_ext,
						   bc_u_ext,
						   ebqe_bc_flux_u_diff_ext[ebNE_kb],
						   mom_u_diff_ten_ext,
						   grad_u_ext,
						   u_ext,
						   h_penalty,//ebqe_penalty_ext[ebNE_kb],
						   flux_mom_u_diff_ext);
	  RANS2PV2::exteriorNumericalDiffusiveFlux(eps_rho,
						   ebqe_phi_ext[ebNE_kb],
						   sdInfo_v_v_rowptr,
						   sdInfo_v_v_colind,
						   isDOFBoundary_v[ebNE_kb],
						   isDiffusiveFluxBoundary_v[ebNE_kb],
						   normal,
						   bc_mom_v_diff_ten_ext,
						   bc_v_ext,
						   ebqe_bc_flux_v_diff_ext[ebNE_kb],
						   mom_v_diff_ten_ext,
						   grad_v_ext,
						   v_ext,
						   h_penalty,//ebqe_penalty_ext[ebNE_kb],
						   flux_mom_v_diff_ext);
	  RANS2PV2::exteriorNumericalDiffusiveFlux(eps_rho,
						   ebqe_phi_ext[ebNE_kb],
						   sdInfo_w_w_rowptr,
						   sdInfo_w_w_colind,
						   isDOFBoundary_w[ebNE_kb],
						   isDiffusiveFluxBoundary_w[ebNE_kb],
						   normal,
						   bc_mom_w_diff_ten_ext,
						   bc_w_ext,
						   ebqe_bc_flux_w_diff_ext[ebNE_kb],
						   mom_w_diff_ten_ext,
						   grad_w_ext,
						   w_ext,
						   h_penalty,//ebqe_penalty_ext[ebNE_kb],
						   flux_mom_w_diff_ext);
	  flux[ebN*nQuadraturePoints_elementBoundary+kb] = flux_mass_ext;
	  //flux[ebN*nQuadraturePoints_elementBoundary+kb] = 0.0;//cek debug
	  //
	  //update residuals
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_p[i] += ck.ExteriorElementBoundaryFlux(flux_mass_ext,p_test_dS[i]);
	      globalConservationError += ck.ExteriorElementBoundaryFlux(flux_mass_ext,p_test_dS[i]);

	      elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_mom_u_adv_ext,vel_test_dS[i])+
	      	ck.ExteriorElementBoundaryFlux(flux_mom_u_diff_ext,vel_test_dS[i]); 

	      elementResidual_v[i] += ck.ExteriorElementBoundaryFlux(flux_mom_v_adv_ext,vel_test_dS[i]) +
	      	ck.ExteriorElementBoundaryFlux(flux_mom_v_diff_ext,vel_test_dS[i]); 
	       
	      elementResidual_w[i] += ck.ExteriorElementBoundaryFlux(flux_mom_w_adv_ext,vel_test_dS[i]) +
	      	ck.ExteriorElementBoundaryFlux(flux_mom_w_diff_ext,vel_test_dS[i]); 
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

	  globalResidual[offset_p+stride_p*p_l2g[eN_i]]+=elementResidual_p[i];
	  globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
	  globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
	  globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
	}//i
      // //
      // //debug
      // //
      // for(int i=0;i<nDOF_test_element;i++) 
      //   { 
      // 	  std::cout<<"ebNE "<<ebNE<<" i "<<i<<std::endl;
      // 	  std::cout<<"r_p"<<elementResidual_p[i]<<std::endl;
      // 	  std::cout<<"r_u"<<elementResidual_u[i]<<std::endl;
      // 	  std::cout<<"r_v"<<elementResidual_v[i]<<std::endl;
      // 	  std::cout<<"r_w"<<elementResidual_w[i]<<std::endl;
      // 	}

    }//ebNE
}

extern "C" void calculateJacobian_RANS2PV2(//element
					   double* mesh_trial_ref,
					   double* mesh_grad_trial_ref,
					   double* mesh_dof,
					   int* mesh_l2g,
					   double* dV_ref,
					   double* p_trial_ref,
					   double* p_grad_trial_ref,
					   double* p_test_ref,
					   double* p_grad_test_ref,
					   double* vel_trial_ref,
					   double* vel_grad_trial_ref,
					   double* vel_test_ref,
					   double* vel_grad_test_ref,
					   //element boundary
					   double* mesh_trial_trace_ref,
					   double* mesh_grad_trial_trace_ref,
					   double* dS_ref,
					   double* p_trial_trace_ref,
					   double* p_grad_trial_trace_ref,
					   double* p_test_trace_ref,
					   double* p_grad_test_trace_ref,
					   double* vel_trial_trace_ref,
					   double* vel_grad_trial_trace_ref,
					   double* vel_test_trace_ref,
					   double* vel_grad_test_trace_ref,					 
					   double* normal_ref,
					   double* boundaryJac_ref,
					   //physics
					   int nElements_global,
					   double alphaBDF,
					   double epsFact_rho,
					   double epsFact_mu,
					   double sigma,
					   double rho_0,
					   double nu_0,
					   double rho_1,
					   double nu_1,
					   double Ct_sge,
					   double Cd_sge,
					   double C_dg,
					   int* p_l2g, 
					   int* vel_l2g,
					   double* p_dof, double* u_dof, double* v_dof, double* w_dof,
					   double* g,
					   double* phi,
					   double* normal_phi,
					   double* kappa_phi,
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
					   double* ebqe_phi_ext,
					   double* ebqe_normal_phi_ext,
					   double* ebqe_kappa_phi_ext,
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
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;
  //
  //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
  //
  for(int eN=0;eN<nElements_global;eN++)
    {
      register double eps_rho,eps_mu;

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
	    eN_k_nSpace = eN_k*nSpace,
	    eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	  //declare local storage
	  register double p=0.0,u=0.0,v=0.0,w=0.0,
	    grad_p[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
	    mom_u_acc=0.0,
	    dmom_u_acc_u=0.0,
	    mom_v_acc=0.0,
	    dmom_v_acc_v=0.0,
	    mom_w_acc=0.0,
	    dmom_w_acc_w=0.0,
	    mass_adv[nSpace],
	    dmass_adv_u[nSpace],
	    dmass_adv_v[nSpace],
	    dmass_adv_w[nSpace],
	    mom_u_adv[nSpace],
	    dmom_u_adv_u[nSpace],
	    dmom_u_adv_v[nSpace],
	    dmom_u_adv_w[nSpace],
	    mom_v_adv[nSpace],
	    dmom_v_adv_u[nSpace],
	    dmom_v_adv_v[nSpace],
	    dmom_v_adv_w[nSpace],
	    mom_w_adv[nSpace],
	    dmom_w_adv_u[nSpace],
	    dmom_w_adv_v[nSpace],
	    dmom_w_adv_w[nSpace],
	    mom_u_diff_ten[nSpace],
	    mom_v_diff_ten[nSpace],
	    mom_w_diff_ten[nSpace],
	    mom_uv_diff_ten[1],
	    mom_uw_diff_ten[1],
	    mom_vu_diff_ten[1],
	    mom_vw_diff_ten[1],
	    mom_wu_diff_ten[1],
	    mom_wv_diff_ten[1],
	    mom_u_source=0.0,
	    mom_v_source=0.0,
	    mom_w_source=0.0,
	    mom_u_ham=0.0,
	    dmom_u_ham_grad_p[nSpace],
	    mom_v_ham=0.0,
	    dmom_v_ham_grad_p[nSpace],
	    mom_w_ham=0.0,
	    dmom_w_ham_grad_p[nSpace],
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
	    dpdeResidual_p_u[nDOF_trial_element],dpdeResidual_p_v[nDOF_trial_element],dpdeResidual_p_w[nDOF_trial_element],
	    dpdeResidual_u_p[nDOF_trial_element],dpdeResidual_u_u[nDOF_trial_element],
	    dpdeResidual_v_p[nDOF_trial_element],dpdeResidual_v_v[nDOF_trial_element],
	    dpdeResidual_w_p[nDOF_trial_element],dpdeResidual_w_w[nDOF_trial_element],
	    Lstar_u_p[nDOF_test_element],
	    Lstar_v_p[nDOF_test_element],
	    Lstar_w_p[nDOF_test_element],
	    Lstar_u_u[nDOF_test_element],
	    Lstar_v_v[nDOF_test_element],
	    Lstar_w_w[nDOF_test_element],
	    Lstar_p_u[nDOF_test_element],
	    Lstar_p_v[nDOF_test_element],
	    Lstar_p_w[nDOF_test_element],
	    subgridError_p=0.0,
	    subgridError_u=0.0,
	    subgridError_v=0.0,
	    subgridError_w=0.0,	    
	    dsubgridError_p_u[nDOF_trial_element],
	    dsubgridError_p_v[nDOF_trial_element],
	    dsubgridError_p_w[nDOF_trial_element],
	    dsubgridError_u_p[nDOF_trial_element],
	    dsubgridError_u_u[nDOF_trial_element],
	    dsubgridError_v_p[nDOF_trial_element],
	    dsubgridError_v_v[nDOF_trial_element],
	    dsubgridError_w_p[nDOF_trial_element],
	    dsubgridError_w_w[nDOF_trial_element],
	    tau_p=0.0,
	    tau_v=0.0,
	    jac[nSpace*nSpace],
	    jacDet,
	    jacInv[nSpace*nSpace],
	    p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
	    dV,
	    p_test_dV[nDOF_test_element],vel_test_dV[nDOF_test_element],
	    p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
	    x,y,z,
	    G[nSpace*nSpace],G_dd_G,tr_G,h_phi, vel[nSpace];
	  //get jacobian, etc for mapping reference element
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
	  ck.calculateG(jacInv,G,G_dd_G,tr_G);
	  ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);
	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;
	  //get the trial function gradients
	  ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
	  ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	  //get the solution 	
	  ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p);
	  ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],u);
	  ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],v);
	  ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w);
	  //get the solution gradients
	  ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p);
	  ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_u);
	  ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_v);
	  ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV;
	      vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
	      for (int I=0;I<nSpace;I++)
		{
		  p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		  vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin}
		}
	    }
	  //cek debug
	  // q_velocity_last[eN_k_nSpace+0]=1.0;
	  // q_velocity_last[eN_k_nSpace+1]=1.0;
	  // q_velocity_last[eN_k_nSpace+2]=1.0;
	  // //
	  // //debugging section for finite element calculations on interior
	  // //
	  // std::cout<<"eN = "<<eN<<" k = "<<k<<std::endl;
	  // for (int j=0;j<nDOF_trial_element;j++)
	  //   {
	  //     std::cout<<"p_trial["<<j<<"] "<<p_trial_ref[k*nDOF_test_element+j]<<std::endl;
	  //     std::cout<<"vel_trial["<<j<<"] "<<vel_trial_ref[k*nDOF_test_element+j]<<std::endl;
	  //     for (int I=0;I<nSpace;I++)
	  // 	{
	  // 	  std::cout<<"p_grad_trial["<<j<<","<<I<<"] "<<p_grad_trial[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"vel_grad_trial["<<j<<","<<I<<"] "<<vel_grad_trial[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"p_grad_test_dV["<<j<<","<<I<<"] "<<p_grad_test_dV[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"vel_grad_test_dV["<<j<<","<<I<<"] "<<vel_grad_test_dV[j*nSpace+I]<<std::endl;
	  // 	}
	  //   }
	  // std::cout<<"p "<<p<<std::endl;
	  // std::cout<<"u "<<u<<std::endl;
	  // std::cout<<"v "<<v<<std::endl;
	  // std::cout<<"w "<<w<<std::endl;
	  // for(int I=0;I<nSpace;I++)
	  //   {
	  //     std::cout<<"grad_p["<<I<<"] "<<grad_p[I]<<std::endl;
	  //     std::cout<<"grad_u["<<I<<"] "<<grad_u[I]<<std::endl;
	  //     std::cout<<"grad_v["<<I<<"] "<<grad_v[I]<<std::endl;
	  //     std::cout<<"grad_w["<<I<<"] "<<grad_w[I]<<std::endl;
	  //   }
          //
          //calculate pde coefficients and derivatives at quadrature points
          //
	  //cek debug
	  //eps_rho=0.1;
	  //eps_mu=0.1;
	  RANS2PV2::evaluateCoefficients(eps_rho,
					 eps_mu,
					 sigma,
					 rho_0,
					 nu_0,
					 rho_1,
					 nu_1,
					 g,
					 phi[eN_k],
					 &normal_phi[eN_k_nSpace],
					 kappa_phi[eN_k],
					 p,
					 grad_p,
					 u,
					 v,
					 w,
					 mom_u_acc,
					 dmom_u_acc_u,
					 mom_v_acc,
					 dmom_v_acc_v,
					 mom_w_acc,
					 dmom_w_acc_w,
					 mass_adv,
					 dmass_adv_u,
					 dmass_adv_v,
					 dmass_adv_w,
					 mom_u_adv,
					 dmom_u_adv_u,
					 dmom_u_adv_v,
					 dmom_u_adv_w,
					 mom_v_adv,
					 dmom_v_adv_u,
					 dmom_v_adv_v,
					 dmom_v_adv_w,
					 mom_w_adv,
					 dmom_w_adv_u,
					 dmom_w_adv_v,
					 dmom_w_adv_w,
					 mom_u_diff_ten,
					 mom_v_diff_ten,
					 mom_w_diff_ten,
					 mom_uv_diff_ten,
					 mom_uw_diff_ten,
					 mom_vu_diff_ten,
					 mom_vw_diff_ten,
					 mom_wu_diff_ten,
					 mom_wv_diff_ten,
					 mom_u_source,
					 mom_v_source,
					 mom_w_source,
					 mom_u_ham,
					 dmom_u_ham_grad_p,
					 mom_v_ham,
					 dmom_v_ham_grad_p,
					 mom_w_ham,
					 dmom_w_ham_grad_p);          
          //
          //moving mesh
          //
          //omit for now
          //
          //calculate time derivatives
          //
          ck.bdf(alphaBDF,
		 q_mom_u_acc_beta_bdf[eN_k],
		 mom_u_acc,
		 dmom_u_acc_u,
		 mom_u_acc_t,
		 dmom_u_acc_u_t);
          ck.bdf(alphaBDF,
		 q_mom_v_acc_beta_bdf[eN_k],
		 mom_v_acc,
		 dmom_v_acc_v,
		 mom_v_acc_t,
		 dmom_v_acc_v_t);
          ck.bdf(alphaBDF,
		 q_mom_w_acc_beta_bdf[eN_k],
		 mom_w_acc,
		 dmom_w_acc_w,
		 mom_w_acc_t,
		 dmom_w_acc_w_t);
          //
          //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
          //
          //calculate the Jacobian of strong residual
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      register int j_nSpace = j*nSpace;
	      dpdeResidual_p_u[j]=ck.AdvectionJacobian_strong(dmass_adv_u,&vel_grad_trial[j_nSpace]);
	      dpdeResidual_p_v[j]=ck.AdvectionJacobian_strong(dmass_adv_v,&vel_grad_trial[j_nSpace]);
	      dpdeResidual_p_w[j]=ck.AdvectionJacobian_strong(dmass_adv_w,&vel_grad_trial[j_nSpace]);

	      dpdeResidual_u_p[j]=ck.HamiltonianJacobian_strong(dmom_u_ham_grad_p,&p_grad_trial[j_nSpace]);
	      dpdeResidual_u_u[j]=ck.MassJacobian_strong(dmom_u_acc_u_t,vel_trial_ref[k*nDOF_trial_element+j]) +
		ck.AdvectionJacobian_strong(&q_velocity_last[eN_k_nSpace],&vel_grad_trial[j_nSpace]);
	      
	      dpdeResidual_v_p[j]=ck.HamiltonianJacobian_strong(dmom_v_ham_grad_p,&p_grad_trial[j_nSpace]);
	      dpdeResidual_v_v[j]=ck.MassJacobian_strong(dmom_v_acc_v_t,vel_trial_ref[k*nDOF_trial_element+j]) +
		ck.AdvectionJacobian_strong(&q_velocity_last[eN_k_nSpace],&vel_grad_trial[j_nSpace]);
	      
	      dpdeResidual_w_p[j]=ck.HamiltonianJacobian_strong(dmom_w_ham_grad_p,&p_grad_trial[j_nSpace]);
	      dpdeResidual_w_w[j]=ck.MassJacobian_strong(dmom_w_acc_w_t,vel_trial_ref[k*nDOF_trial_element+j]) + 
		ck.AdvectionJacobian_strong(&q_velocity_last[eN_k_nSpace],&vel_grad_trial[j_nSpace]);
            }
          //calculate tau and tau*Res
	  RANS2PV2::calculateSubgridError_tau(Ct_sge,
					      Cd_sge,
					      G,
					      G_dd_G,
					      tr_G,
					      dmom_u_acc_u,//rho
					      dmom_u_acc_u_t/dmom_u_acc_u,//Dt
	  				      &q_velocity_last[eN_k_nSpace],//v
					      mom_u_diff_ten[1],//mu
	  				      tau_p,
					      tau_v,
					      q_cfl[eN_k]);
	  //cek debug
	  //tau_p = 0.0;
	  //tau_v = 0.0;
	  RANS2PV2::calculateSubgridErrorDerivatives_tauRes(tau_p,
							    tau_v,
							    dpdeResidual_p_u,
							    dpdeResidual_p_v,
							    dpdeResidual_p_w,
							    dpdeResidual_u_p,
							    dpdeResidual_u_u,
							    dpdeResidual_v_p,
							    dpdeResidual_v_v,
							    dpdeResidual_w_p,
							    dpdeResidual_w_w,
							    dsubgridError_p_u,
							    dsubgridError_p_v,
							    dsubgridError_p_w,
							    dsubgridError_u_p,
							    dsubgridError_u_u,
							    dsubgridError_v_p,
							    dsubgridError_v_v,
							    dsubgridError_w_p,
							    dsubgridError_w_w);
	  //cek/ido todo add options for lagging velocity and turning off RBLES terms
         // Velocity  ==> SUPG   
	    //vel[0] = u;
	    //vel[1] = v;
	    //vel[2] = w;
 
        // Velocity + subgrid ==> RBLES
	    vel[0] = u + subgridError_u;
	    vel[1] = v + subgridError_v;
	    vel[2] = w + subgridError_w;

          //calculate the adjoint times the test functions
          for (int i=0;i<nDOF_test_element;i++)
            {
	      register int i_nSpace = i*nSpace;
	      Lstar_u_p[i]=ck.Advection_adjoint(dmass_adv_u,&p_grad_test_dV[i_nSpace]);
	      Lstar_v_p[i]=ck.Advection_adjoint(dmass_adv_v,&p_grad_test_dV[i_nSpace]);
	      Lstar_w_p[i]=ck.Advection_adjoint(dmass_adv_w,&p_grad_test_dV[i_nSpace]);
	      Lstar_u_u[i]=ck.Advection_adjoint(vel,&vel_grad_test_dV[i_nSpace]);
	      Lstar_v_v[i]=ck.Advection_adjoint(vel,&vel_grad_test_dV[i_nSpace]);
	      Lstar_w_w[i]=ck.Advection_adjoint(vel,&vel_grad_test_dV[i_nSpace]);
	      Lstar_p_u[i]=ck.Hamiltonian_adjoint(dmom_u_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
	      Lstar_p_v[i]=ck.Hamiltonian_adjoint(dmom_v_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
	      Lstar_p_w[i]=ck.Hamiltonian_adjoint(dmom_w_ham_grad_p,&vel_grad_test_dV[i_nSpace]);
            }


	  //cek todo add RBLES terms consistent to residual modifications or ignore them partials w.r.t the additional RBLES terms
  	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      register int i_nSpace = i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  register int j_nSpace = j*nSpace;
		  elementJacobian_p_p[i][j] += ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_p[i]) + 
		    ck.SubgridErrorJacobian(dsubgridError_v_p[j],Lstar_v_p[i]) + 
		    ck.SubgridErrorJacobian(dsubgridError_w_p[j],Lstar_w_p[i]); 

		  elementJacobian_p_u[i][j] += ck.AdvectionJacobian_weak(dmass_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u_p[i]); 
		  elementJacobian_p_v[i][j] += ck.AdvectionJacobian_weak(dmass_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_v_v[j],Lstar_v_p[i]); 
		  elementJacobian_p_w[i][j] += ck.AdvectionJacobian_weak(dmass_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_w_w[j],Lstar_w_p[i]); 

		  elementJacobian_u_p[i][j] += ck.HamiltonianJacobian_weak(dmom_u_ham_grad_p,&p_grad_trial[j_nSpace],vel_test_dV[i]) + 
		    ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_u[i]); 
		  elementJacobian_u_u[i][j] += ck.MassJacobian_weak(dmom_u_acc_u_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
		    ck.AdvectionJacobian_weak(dmom_u_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
		    ck.SimpleDiffusionJacobian_weak(sdInfo_u_u_rowptr,sdInfo_u_u_colind,mom_u_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_u[i]) + 
		    ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u_u[i]) + 
		    ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]); 
		  elementJacobian_u_v[i][j] += ck.AdvectionJacobian_weak(dmom_u_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SimpleDiffusionJacobian_weak(sdInfo_u_v_rowptr,sdInfo_u_v_colind,mom_uv_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_u[i]); 
		  elementJacobian_u_w[i][j] += ck.AdvectionJacobian_weak(dmom_u_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SimpleDiffusionJacobian_weak(sdInfo_u_w_rowptr,sdInfo_u_w_colind,mom_uw_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_w[j],Lstar_p_u[i]); 

		  elementJacobian_v_p[i][j] += ck.HamiltonianJacobian_weak(dmom_v_ham_grad_p,&p_grad_trial[j_nSpace],vel_test_dV[i]) + 
		    ck.SubgridErrorJacobian(dsubgridError_v_p[j],Lstar_v_v[i]); 
		  elementJacobian_v_u[i][j] += ck.AdvectionJacobian_weak(dmom_v_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SimpleDiffusionJacobian_weak(sdInfo_v_u_rowptr,sdInfo_v_u_colind,mom_vu_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_v[i]);
		  elementJacobian_v_v[i][j] += ck.MassJacobian_weak(dmom_v_acc_v_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
		    ck.AdvectionJacobian_weak(dmom_v_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
		    ck.SimpleDiffusionJacobian_weak(sdInfo_v_v_rowptr,sdInfo_v_v_colind,mom_v_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_v[i]) +
		    ck.SubgridErrorJacobian(dsubgridError_v_v[j],Lstar_v_v[i]) + 
		    ck.NumericalDiffusionJacobian(q_numDiff_v_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]); 
		  elementJacobian_v_w[i][j] += ck.AdvectionJacobian_weak(dmom_v_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +  
		    ck.SimpleDiffusionJacobian_weak(sdInfo_v_w_rowptr,sdInfo_v_w_colind,mom_vw_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_w[j],Lstar_p_v[i]);

		  elementJacobian_w_p[i][j] += ck.HamiltonianJacobian_weak(dmom_w_ham_grad_p,&p_grad_trial[j_nSpace],vel_test_dV[i]) + 
		    ck.SubgridErrorJacobian(dsubgridError_w_p[j],Lstar_w_w[i]); 
		  elementJacobian_w_u[i][j] += ck.AdvectionJacobian_weak(dmom_w_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +  
		    ck.SimpleDiffusionJacobian_weak(sdInfo_w_u_rowptr,sdInfo_w_u_colind,mom_wu_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_w[i]); 
		  elementJacobian_w_v[i][j] += ck.AdvectionJacobian_weak(dmom_w_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SimpleDiffusionJacobian_weak(sdInfo_w_v_rowptr,sdInfo_w_v_colind,mom_wv_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_w[i]); 
		  elementJacobian_w_w[i][j] += ck.MassJacobian_weak(dmom_w_acc_w_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
		    ck.AdvectionJacobian_weak(dmom_w_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +  
		    ck.SimpleDiffusionJacobian_weak(sdInfo_w_w_rowptr,sdInfo_w_w_colind,mom_w_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
		    ck.SubgridErrorJacobian(dsubgridError_p_w[j],Lstar_p_w[i]) + 
		    ck.SubgridErrorJacobian(dsubgridError_w_w[j],Lstar_w_w[i]) + 
		    ck.NumericalDiffusionJacobian(q_numDiff_w_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]); 
		}//j
            }//i
	}//k
      //
      //load into element Jacobian into global Jacobian
      //
      for (int i=0;i<nDOF_test_element;i++)
	{
	  register int eN_i = eN*nDOF_test_element+i;
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      register int eN_i_j = eN_i*nDOF_trial_element+j;
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
      // //
      // //debug element jacobian
      // //
      // std::cout<<"element jacobian"<<std::endl;
      // for (int i=0;i<nDOF_test_element;i++)
      // 	{
      // 	  for (int j=0;j<nDOF_trial_element;j++)
      // 	    {
      // 	      std::cout << elementJacobian_p_p[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_p_u[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_p_v[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_p_w[i][j]<<std::endl;

      // 	      std::cout << elementJacobian_u_p[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_u_u[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_u_v[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_u_w[i][j]<<std::endl;

      // 	      std::cout << elementJacobian_v_p[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_v_u[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_v_v[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_v_w[i][j]<<std::endl;

      // 	      std::cout << elementJacobian_w_p[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_w_u[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_w_v[i][j]<<std::endl;
      // 	      std::cout << elementJacobian_w_w[i][j]<<std::endl;
      // 	    }//j
      // 	}//i
    }//elements
  //
  //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
  //
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      register int ebN = exteriorElementBoundariesArray[ebNE],
	eN  = elementBoundaryElementsArray[ebN*2+0],
	eN_nDOF_trial_element = eN*nDOF_trial_element,
	ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
      register double eps_rho,eps_mu;
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace,
	    ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
	    ebN_local_kb_nSpace = ebN_local_kb*nSpace;

	  register double p_ext=0.0,
	    u_ext=0.0,
	    v_ext=0.0,
	    w_ext=0.0,
	    grad_p_ext[nSpace],
	    grad_u_ext[nSpace],
	    grad_v_ext[nSpace],
	    grad_w_ext[nSpace],
	    mom_u_acc_ext=0.0,
	    dmom_u_acc_u_ext=0.0,
	    mom_v_acc_ext=0.0,
	    dmom_v_acc_v_ext=0.0,
	    mom_w_acc_ext=0.0,
	    dmom_w_acc_w_ext=0.0,
	    mass_adv_ext[nSpace],
	    dmass_adv_u_ext[nSpace],
	    dmass_adv_v_ext[nSpace],
	    dmass_adv_w_ext[nSpace],
	    mom_u_adv_ext[nSpace],
	    dmom_u_adv_u_ext[nSpace],
	    dmom_u_adv_v_ext[nSpace],
	    dmom_u_adv_w_ext[nSpace],
	    mom_v_adv_ext[nSpace],
	    dmom_v_adv_u_ext[nSpace],
	    dmom_v_adv_v_ext[nSpace],
	    dmom_v_adv_w_ext[nSpace],
	    mom_w_adv_ext[nSpace],
	    dmom_w_adv_u_ext[nSpace],
	    dmom_w_adv_v_ext[nSpace],
	    dmom_w_adv_w_ext[nSpace],
	    mom_u_diff_ten_ext[nSpace],
	    mom_v_diff_ten_ext[nSpace],
	    mom_w_diff_ten_ext[nSpace],
	    mom_uv_diff_ten_ext[1],
	    mom_uw_diff_ten_ext[1],
	    mom_vu_diff_ten_ext[1],
	    mom_vw_diff_ten_ext[1],
	    mom_wu_diff_ten_ext[1],
	    mom_wv_diff_ten_ext[1],
	    mom_u_source_ext=0.0,
	    mom_v_source_ext=0.0,
	    mom_w_source_ext=0.0,
	    mom_u_ham_ext=0.0,
	    dmom_u_ham_grad_p_ext[nSpace],
	    mom_v_ham_ext=0.0,
	    dmom_v_ham_grad_p_ext[nSpace],
	    mom_w_ham_ext=0.0,
	    dmom_w_ham_grad_p_ext[nSpace],
	    dmom_u_adv_p_ext[nSpace],
	    dmom_v_adv_p_ext[nSpace],
	    dmom_w_adv_p_ext[nSpace],
	    dflux_mass_u_ext=0.0,
	    dflux_mass_v_ext=0.0,
	    dflux_mass_w_ext=0.0,
	    dflux_mom_u_adv_p_ext=0.0,
	    dflux_mom_u_adv_u_ext=0.0,
	    dflux_mom_u_adv_v_ext=0.0,
	    dflux_mom_u_adv_w_ext=0.0,
	    dflux_mom_v_adv_p_ext=0.0,
	    dflux_mom_v_adv_u_ext=0.0,
	    dflux_mom_v_adv_v_ext=0.0,
	    dflux_mom_v_adv_w_ext=0.0,
	    dflux_mom_w_adv_p_ext=0.0,
	    dflux_mom_w_adv_u_ext=0.0,
	    dflux_mom_w_adv_v_ext=0.0,
	    dflux_mom_w_adv_w_ext=0.0,
	    bc_p_ext=0.0,
	    bc_u_ext=0.0,
	    bc_v_ext=0.0,
	    bc_w_ext=0.0,
	    bc_mom_u_acc_ext=0.0,
	    bc_dmom_u_acc_u_ext=0.0,
	    bc_mom_v_acc_ext=0.0,
	    bc_dmom_v_acc_v_ext=0.0,
	    bc_mom_w_acc_ext=0.0,
	    bc_dmom_w_acc_w_ext=0.0,
	    bc_mass_adv_ext[nSpace],
	    bc_dmass_adv_u_ext[nSpace],
	    bc_dmass_adv_v_ext[nSpace],
	    bc_dmass_adv_w_ext[nSpace],
	    bc_mom_u_adv_ext[nSpace],
	    bc_dmom_u_adv_u_ext[nSpace],
	    bc_dmom_u_adv_v_ext[nSpace],
	    bc_dmom_u_adv_w_ext[nSpace],
	    bc_mom_v_adv_ext[nSpace],
	    bc_dmom_v_adv_u_ext[nSpace],
	    bc_dmom_v_adv_v_ext[nSpace],
	    bc_dmom_v_adv_w_ext[nSpace],
	    bc_mom_w_adv_ext[nSpace],
	    bc_dmom_w_adv_u_ext[nSpace],
	    bc_dmom_w_adv_v_ext[nSpace],
	    bc_dmom_w_adv_w_ext[nSpace],
	    bc_mom_u_diff_ten_ext[nSpace],
	    bc_mom_v_diff_ten_ext[nSpace],
	    bc_mom_w_diff_ten_ext[nSpace],
	    bc_mom_uv_diff_ten_ext[1],
	    bc_mom_uw_diff_ten_ext[1],
	    bc_mom_vu_diff_ten_ext[1],
	    bc_mom_vw_diff_ten_ext[1],
	    bc_mom_wu_diff_ten_ext[1],
	    bc_mom_wv_diff_ten_ext[1],
	    bc_mom_u_source_ext=0.0,
	    bc_mom_v_source_ext=0.0,
	    bc_mom_w_source_ext=0.0,
	    bc_mom_u_ham_ext=0.0,
	    bc_dmom_u_ham_grad_p_ext[nSpace],
	    bc_mom_v_ham_ext=0.0,
	    bc_dmom_v_ham_grad_p_ext[nSpace],
	    bc_mom_w_ham_ext=0.0,
	    bc_dmom_w_ham_grad_p_ext[nSpace],
	    fluxJacobian_p_p[nDOF_trial_element],
	    fluxJacobian_p_u[nDOF_trial_element],
	    fluxJacobian_p_v[nDOF_trial_element],
	    fluxJacobian_p_w[nDOF_trial_element],
	    fluxJacobian_u_p[nDOF_trial_element],
	    fluxJacobian_u_u[nDOF_trial_element],
	    fluxJacobian_u_v[nDOF_trial_element],
	    fluxJacobian_u_w[nDOF_trial_element],
	    fluxJacobian_v_p[nDOF_trial_element],
	    fluxJacobian_v_u[nDOF_trial_element],
	    fluxJacobian_v_v[nDOF_trial_element],
	    fluxJacobian_v_w[nDOF_trial_element],
	    fluxJacobian_w_p[nDOF_trial_element],
	    fluxJacobian_w_u[nDOF_trial_element],
	    fluxJacobian_w_v[nDOF_trial_element],
	    fluxJacobian_w_w[nDOF_trial_element],
	    jac_ext[nSpace*nSpace],
	    jacDet_ext,
	    jacInv_ext[nSpace*nSpace],
	    boundaryJac[nSpace*(nSpace-1)],
	    metricTensor[(nSpace-1)*(nSpace-1)],
	    metricTensorDetSqrt,
	    p_grad_trial_trace[nDOF_trial_element*nSpace],
	    vel_grad_trial_trace[nDOF_trial_element*nSpace],
	    dS,
	    p_test_dS[nDOF_test_element],
	    vel_test_dS[nDOF_test_element],
	    normal[3],
	    x_ext,y_ext,z_ext,
	    G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty;
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
	  dS = metricTensorDetSqrt*dS_ref[kb];
	  ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	  ck.calculateGScale(G,&ebqe_normal_phi_ext[ebNE_kb_nSpace],h_phi);
	  eps_rho = epsFact_rho*h_phi;
	  eps_mu = epsFact_mu*h_phi;
	  //compute shape and solution information
	  //shape
	  ck.gradTrialFromRef(&p_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
	  ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace);
	  //solution and gradients	
	  ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],p_ext);
	  ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	  ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext);
	  ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],w_ext);
	  ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
	  ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_u_ext);
	  ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_v_ext);
	  ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_w_ext);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      p_test_dS[j] = p_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	      vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
	    }
	  // //
	  // //debugging section for finite element calculations on exterior
	  // //
	  // std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl;
	  // for (int j=0;j<nDOF_trial_element;j++)
	  //   {
	  //     std::cout<<"p_trial_trace["<<j<<"] "<<p_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	  //     std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	  //     std::cout<<"p_test_dS["<<j<<"] "<<p_test_dS[j]<<std::endl;
	  //     std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl;
	  //     for (int I=0;I<nSpace;I++)
	  // 	{
	  // 	  std::cout<<"p_grad_trial_trace["<<j<<","<<I<<"] "<<p_grad_trial_trace[j*nSpace+I]<<std::endl;
	  // 	  std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl;
	  // 	}
	  //   }
	  // std::cout<<"p_ext "<<p_ext<<std::endl;
	  // std::cout<<"u_ext "<<u_ext<<std::endl;
	  // std::cout<<"v_ext "<<v_ext<<std::endl;
	  // std::cout<<"w_ext "<<w_ext<<std::endl;
	  // for(int I=0;I<nSpace;I++)
	  //   {
	  //     std::cout<<"grad_p_ext["<<I<<"] "<<grad_p_ext[I]<<std::endl;
	  //     std::cout<<"grad_u_ext["<<I<<"] "<<grad_u_ext[I]<<std::endl;
	  //     std::cout<<"grad_v_ext["<<I<<"] "<<grad_v_ext[I]<<std::endl;
	  //     std::cout<<"grad_w_ext["<<I<<"] "<<grad_w_ext[I]<<std::endl;
	  //   }
	  //
	  //load the boundary values
	  //
	  bc_p_ext = isDOFBoundary_p[ebNE_kb]*ebqe_bc_p_ext[ebNE_kb]+(1-isDOFBoundary_p[ebNE_kb])*p_ext;
	  bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	  bc_v_ext = isDOFBoundary_v[ebNE_kb]*ebqe_bc_v_ext[ebNE_kb]+(1-isDOFBoundary_v[ebNE_kb])*v_ext;
	  bc_w_ext = isDOFBoundary_w[ebNE_kb]*ebqe_bc_w_ext[ebNE_kb]+(1-isDOFBoundary_w[ebNE_kb])*w_ext;
	  // 
	  //calculate the internal and external trace of the pde coefficients 
	  // 
	  //cek debug
	  //eps_rho=0.1;
	  //eps_mu=0.1;
	  RANS2PV2::evaluateCoefficients(eps_rho,
					 eps_mu,
					 sigma,
					 rho_0,
					 nu_0,
					 rho_1,
					 nu_1,
					 g,
					 ebqe_phi_ext[ebNE_kb],
					 &ebqe_normal_phi_ext[ebNE_kb_nSpace],
					 ebqe_kappa_phi_ext[ebNE_kb],
					 p_ext,
					 grad_p_ext,
					 u_ext,
					 v_ext,
					 w_ext,
					 mom_u_acc_ext,
					 dmom_u_acc_u_ext,
					 mom_v_acc_ext,
					 dmom_v_acc_v_ext,
					 mom_w_acc_ext,
					 dmom_w_acc_w_ext,
					 mass_adv_ext,
					 dmass_adv_u_ext,
					 dmass_adv_v_ext,
					 dmass_adv_w_ext,
					 mom_u_adv_ext,
					 dmom_u_adv_u_ext,
					 dmom_u_adv_v_ext,
					 dmom_u_adv_w_ext,
					 mom_v_adv_ext,
					 dmom_v_adv_u_ext,
					 dmom_v_adv_v_ext,
					 dmom_v_adv_w_ext,
					 mom_w_adv_ext,
					 dmom_w_adv_u_ext,
					 dmom_w_adv_v_ext,
					 dmom_w_adv_w_ext,
					 mom_u_diff_ten_ext,
					 mom_v_diff_ten_ext,
					 mom_w_diff_ten_ext,
					 mom_uv_diff_ten_ext,
					 mom_uw_diff_ten_ext,
					 mom_vu_diff_ten_ext,
					 mom_vw_diff_ten_ext,
					 mom_wu_diff_ten_ext,
					 mom_wv_diff_ten_ext,
					 mom_u_source_ext,
					 mom_v_source_ext,
					 mom_w_source_ext,
					 mom_u_ham_ext,
					 dmom_u_ham_grad_p_ext,
					 mom_v_ham_ext,
					 dmom_v_ham_grad_p_ext,
					 mom_w_ham_ext,
					 dmom_w_ham_grad_p_ext);          
	  RANS2PV2::evaluateCoefficients(eps_rho,
					 eps_mu,
					 sigma,
					 rho_0,
					 nu_0,
					 rho_1,
					 nu_1,
					 g,
					 ebqe_phi_ext[ebNE_kb],
					 &ebqe_normal_phi_ext[ebNE_kb_nSpace],
					 ebqe_kappa_phi_ext[ebNE_kb],
					 bc_p_ext,
					 grad_p_ext, //cek shouldn't be used
					 bc_u_ext,
					 bc_v_ext,
					 bc_w_ext,
					 bc_mom_u_acc_ext,
					 bc_dmom_u_acc_u_ext,
					 bc_mom_v_acc_ext,
					 bc_dmom_v_acc_v_ext,
					 bc_mom_w_acc_ext,
					 bc_dmom_w_acc_w_ext,
					 bc_mass_adv_ext,
					 bc_dmass_adv_u_ext,
					 bc_dmass_adv_v_ext,
					 bc_dmass_adv_w_ext,
					 bc_mom_u_adv_ext,
					 bc_dmom_u_adv_u_ext,
					 bc_dmom_u_adv_v_ext,
					 bc_dmom_u_adv_w_ext,
					 bc_mom_v_adv_ext,
					 bc_dmom_v_adv_u_ext,
					 bc_dmom_v_adv_v_ext,
					 bc_dmom_v_adv_w_ext,
					 bc_mom_w_adv_ext,
					 bc_dmom_w_adv_u_ext,
					 bc_dmom_w_adv_v_ext,
					 bc_dmom_w_adv_w_ext,
					 bc_mom_u_diff_ten_ext,
					 bc_mom_v_diff_ten_ext,
					 bc_mom_w_diff_ten_ext,
					 bc_mom_uv_diff_ten_ext,
					 bc_mom_uw_diff_ten_ext,
					 bc_mom_vu_diff_ten_ext,
					 bc_mom_vw_diff_ten_ext,
					 bc_mom_wu_diff_ten_ext,
					 bc_mom_wv_diff_ten_ext,
					 bc_mom_u_source_ext,
					 bc_mom_v_source_ext,
					 bc_mom_w_source_ext,
					 bc_mom_u_ham_ext,
					 bc_dmom_u_ham_grad_p_ext,
					 bc_mom_v_ham_ext,
					 bc_dmom_v_ham_grad_p_ext,
					 bc_mom_w_ham_ext,
					 bc_dmom_w_ham_grad_p_ext);          
	  // 
	  //calculate the numerical fluxes 
	  // 
	  RANS2PV2::exteriorNumericalAdvectiveFluxDerivatives(isDOFBoundary_p[ebNE_kb],
							      isDOFBoundary_u[ebNE_kb],
							      isDOFBoundary_v[ebNE_kb],
							      isDOFBoundary_w[ebNE_kb],
							      isAdvectiveFluxBoundary_p[ebNE_kb],
							      isAdvectiveFluxBoundary_u[ebNE_kb],
							      isAdvectiveFluxBoundary_v[ebNE_kb],
							      isAdvectiveFluxBoundary_w[ebNE_kb],
							      normal,
							      bc_p_ext,
							      bc_mass_adv_ext,
							      bc_mom_u_adv_ext,
							      bc_mom_v_adv_ext,
							      bc_mom_w_adv_ext,
							      ebqe_bc_flux_mass_ext[ebNE_kb],
							      ebqe_bc_flux_mom_u_adv_ext[ebNE_kb],
							      ebqe_bc_flux_mom_v_adv_ext[ebNE_kb],
							      ebqe_bc_flux_mom_w_adv_ext[ebNE_kb],
							      p_ext,
							      mass_adv_ext,
							      mom_u_adv_ext,
							      mom_v_adv_ext,
							      mom_w_adv_ext,
							      dmass_adv_u_ext,
							      dmass_adv_v_ext,
							      dmass_adv_w_ext,
							      dmom_u_adv_p_ext,
							      dmom_u_adv_u_ext,
							      dmom_u_adv_v_ext,
							      dmom_u_adv_w_ext,
							      dmom_v_adv_p_ext,
							      dmom_v_adv_u_ext,
							      dmom_v_adv_v_ext,
							      dmom_v_adv_w_ext,
							      dmom_w_adv_p_ext,
							      dmom_w_adv_u_ext,
							      dmom_w_adv_v_ext,
							      dmom_w_adv_w_ext,
							      dflux_mass_u_ext,
							      dflux_mass_v_ext,
							      dflux_mass_w_ext,
							      dflux_mom_u_adv_p_ext,
							      dflux_mom_u_adv_u_ext,
							      dflux_mom_u_adv_v_ext,
							      dflux_mom_u_adv_w_ext,
							      dflux_mom_v_adv_p_ext,
							      dflux_mom_v_adv_u_ext,
							      dflux_mom_v_adv_v_ext,
							      dflux_mom_v_adv_w_ext,
							      dflux_mom_w_adv_p_ext,
							      dflux_mom_w_adv_u_ext,
							      dflux_mom_w_adv_v_ext,
							      dflux_mom_w_adv_w_ext);
	  //
	  //calculate the flux jacobian
	  //
	  ck.calculateGScale(G,normal,h_penalty);
	  h_penalty = 10.0/h_penalty;
	  //cek debug
	  //h_penalty = 10.0/0.1;//he = 0.1
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      register int j_nSpace = j*nSpace,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
	      //cek debug
	      //ebqe_penalty_ext[ebNE_kb] = 10.0;
	      //
	      //cek todo add full stress on boundaries

	      fluxJacobian_p_p[j]=0.0;
	      fluxJacobian_p_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_p_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_p_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

	      fluxJacobian_u_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		RANS2PV2::ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
								 ebqe_phi_ext[ebNE_kb],
								 sdInfo_u_u_rowptr,
								 sdInfo_u_u_colind,
								 isDOFBoundary_u[ebNE_kb],
								 normal,
								 mom_u_diff_ten_ext,
								 vel_trial_trace_ref[ebN_local_kb_j],
								 &vel_grad_trial_trace[j_nSpace],
								 h_penalty);//ebqe_penalty_ext[ebNE_kb]);
	      fluxJacobian_u_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_u_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

	      fluxJacobian_v_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_v_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_v_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		RANS2PV2::ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
								 ebqe_phi_ext[ebNE_kb],
								 sdInfo_v_v_rowptr,
								 sdInfo_v_v_colind,
								 isDOFBoundary_v[ebNE_kb],
								 normal,
								 mom_v_diff_ten_ext,
								 vel_trial_trace_ref[ebN_local_kb_j],
								 &vel_grad_trial_trace[j_nSpace],
						 		 h_penalty);//ebqe_penalty_ext[ebNE_kb]);
	      fluxJacobian_v_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

	      fluxJacobian_w_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_w_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_w_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_w_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		RANS2PV2::ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
								 ebqe_phi_ext[ebNE_kb],
								 sdInfo_w_w_rowptr,
								 sdInfo_w_w_colind,
								 isDOFBoundary_w[ebNE_kb],
								 normal,
								 mom_w_diff_ten_ext,
								 vel_trial_trace_ref[ebN_local_kb_j],
								 &vel_grad_trial_trace[j_nSpace],
								 h_penalty);//ebqe_penalty_ext[ebNE_kb]);

	      fluxJacobian_p_p[j]=0.0;
	      fluxJacobian_p_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_p_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_p_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

	      fluxJacobian_u_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		RANS2PV2::ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
								 ebqe_phi_ext[ebNE_kb],
								 sdInfo_u_u_rowptr,
								 sdInfo_u_u_colind,
								 isDOFBoundary_u[ebNE_kb],
								 normal,
								 mom_u_diff_ten_ext,
								 vel_trial_trace_ref[ebN_local_kb_j],
								 &vel_grad_trial_trace[j_nSpace],
								 h_penalty);//ebqe_penalty_ext[ebNE_kb]);
	      fluxJacobian_u_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_u_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

	      fluxJacobian_v_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_v_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_v_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		RANS2PV2::ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
								 ebqe_phi_ext[ebNE_kb],
								 sdInfo_v_v_rowptr,
								 sdInfo_v_v_colind,
								 isDOFBoundary_v[ebNE_kb],
								 normal,
								 mom_v_diff_ten_ext,
								 vel_trial_trace_ref[ebN_local_kb_j],
								 &vel_grad_trial_trace[j_nSpace],
						 		 h_penalty);//ebqe_penalty_ext[ebNE_kb]);
	      fluxJacobian_v_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

	      fluxJacobian_w_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_w_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_w_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
	      fluxJacobian_w_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		RANS2PV2::ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
								 ebqe_phi_ext[ebNE_kb],
								 sdInfo_w_w_rowptr,
								 sdInfo_w_w_colind,
								 isDOFBoundary_w[ebNE_kb],
								 normal,
								 mom_w_diff_ten_ext,
								 vel_trial_trace_ref[ebN_local_kb_j],
								 &vel_grad_trial_trace[j_nSpace],
								 h_penalty);//ebqe_penalty_ext[ebNE_kb]);
	      // //cek debug
	      // fluxJacobian_p_p[j]=0.0;
	      // fluxJacobian_p_u[j]=0.0;
	      // fluxJacobian_p_v[j]=0.0;
	      // fluxJacobian_p_w[j]=0.0;

	      // fluxJacobian_u_p[j]=0.0;
	      // fluxJacobian_u_u[j]=0.0;
	      // fluxJacobian_u_v[j]=0.0;
	      // fluxJacobian_u_w[j]=0.0;

	      // fluxJacobian_v_p[j]=0.0;
	      // fluxJacobian_v_u[j]=0.0;
	      // fluxJacobian_v_v[j]=0.0;
	      // fluxJacobian_v_w[j]=0.0;

	      // fluxJacobian_w_p[j]=0.0;
	      // fluxJacobian_w_u[j]=0.0;
	      // fluxJacobian_w_v[j]=0.0;
	      // fluxJacobian_w_w[j]=0.0;
	      // //cek debug
	    }//j
	  //
	  //update the global Jacobian from the flux Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i = eN*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
		  
		  globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_eb_p_p[ebN_i_j]] += fluxJacobian_p_p[j]*p_test_dS[i];
		  globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_eb_p_u[ebN_i_j]] += fluxJacobian_p_u[j]*p_test_dS[i];
		  globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_eb_p_v[ebN_i_j]] += fluxJacobian_p_v[j]*p_test_dS[i];
		  globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_eb_p_w[ebN_i_j]] += fluxJacobian_p_w[j]*p_test_dS[i];
		   
		  globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_eb_u_p[ebN_i_j]] += fluxJacobian_u_p[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_eb_u_v[ebN_i_j]] += fluxJacobian_u_v[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_eb_u_w[ebN_i_j]] += fluxJacobian_u_w[j]*vel_test_dS[i];
		   
		  globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_eb_v_p[ebN_i_j]] += fluxJacobian_v_p[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_eb_v_u[ebN_i_j]] += fluxJacobian_v_u[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_eb_v_v[ebN_i_j]] += fluxJacobian_v_v[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_eb_v_w[ebN_i_j]] += fluxJacobian_v_w[j]*vel_test_dS[i];
		   
		  globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_eb_w_p[ebN_i_j]] += fluxJacobian_w_p[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_eb_w_u[ebN_i_j]] += fluxJacobian_w_u[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_eb_w_v[ebN_i_j]] += fluxJacobian_w_v[j]*vel_test_dS[i];
		  globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_eb_w_w[ebN_i_j]] += fluxJacobian_w_w[j]*vel_test_dS[i];
		}//j
	    }//i
	  // //debug
	  // std::cout<<"flux jacobian ebNE "<<ebNE<<" kb "<<kb<<std::endl;
	  // for (int i=0;i<nDOF_test_element;i++)
	  //   {
	  //     for (int j=0;j<nDOF_trial_element;j++)
	  // 	{
	  // 	  std::cout<< fluxJacobian_p_p[j]*p_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_p_u[j]*p_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_p_v[j]*p_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_p_w[j]*p_test_dS[i]<<std::endl;
		  
	  // 	  std::cout<< fluxJacobian_u_p[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_u_u[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_u_v[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_u_w[j]*vel_test_dS[i]<<std::endl;
		  
	  // 	  std::cout<< fluxJacobian_v_p[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_v_u[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_v_v[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_v_w[j]*vel_test_dS[i]<<std::endl;
		  
	  // 	  std::cout<< fluxJacobian_w_p[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_w_u[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_w_v[j]*vel_test_dS[i]<<std::endl;
	  // 	  std::cout<< fluxJacobian_w_w[j]*vel_test_dS[i]<<std::endl;
	  // 	}//j
	  //   }//i
	}//kb
    }//ebNE
}//computeJacobian

extern "C" void calculateVelocityAverage_RANS2PV2(int* permutations,
						  int nExteriorElementBoundaries_global,
						  int* exteriorElementBoundariesArray,
						  int nInteriorElementBoundaries_global,
						  int* interiorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  int* elementBoundaryLocalElementBoundariesArray,
						  int* vel_l2g, 
						  double* u_dof, double* v_dof, double* w_dof,
						  double* vel_trial_ref,
						  double* ebqe_velocity,
						  double* velocityAverage)
{
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
    { 
      register int ebN = exteriorElementBoundariesArray[ebNE];
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
	  register int left_kb = permutations[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
					      left_ebN_element*nQuadraturePoints_elementBoundary+
					      kb],
	    right_kb = permutations[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				    right_ebN_element*nQuadraturePoints_elementBoundary+
				    kb];
	  // 
	  //calculate the velocity solution at quadrature points on left and right
	  // 
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      int left_eN_j = left_eN_global*nDOF_trial_element+j;
	      int left_ebN_kb_j = left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		left_kb*nDOF_trial_element +
		j;
	      int right_eN_j = right_eN_global*nDOF_trial_element+j;
	      int right_ebN_kb_j = right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
		right_kb*nDOF_trial_element +
		j;
	      u_left += u_dof[vel_l2g[left_eN_j]]*vel_trial_ref[left_ebN_kb_j]; 
	      v_left += v_dof[vel_l2g[left_eN_j]]*vel_trial_ref[left_ebN_kb_j]; 
	      w_left += w_dof[vel_l2g[left_eN_j]]*vel_trial_ref[left_ebN_kb_j]; 
	      u_right += u_dof[vel_l2g[right_eN_j]]*vel_trial_ref[right_ebN_kb_j]; 
	      v_right += v_dof[vel_l2g[right_eN_j]]*vel_trial_ref[right_ebN_kb_j]; 
	      w_right += w_dof[vel_l2g[right_eN_j]]*vel_trial_ref[right_ebN_kb_j]; 
	    }
	  velocityAverage[ebN_kb_nSpace+0]=0.5*(u_left + u_right);
	  velocityAverage[ebN_kb_nSpace+1]=0.5*(v_left + v_right);
	  velocityAverage[ebN_kb_nSpace+2]=0.5*(w_left + w_right);
	}//ebNI
    }
}
