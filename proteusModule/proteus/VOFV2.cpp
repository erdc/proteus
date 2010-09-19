#include "VOFV2.h"
#include "CompKernel.h"
#include <iostream>

//extern "C" void calculateResidual_VOFV2(int nElements_global,
extern "C" void VOFV2_RES(//element
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
                          int nElements_global,
			double alphaBDF,
			int lag_shockCapturing, /*mwf not used yet*/
			double shockCapturingDiffusion,
			int* u_l2g, 
			double* elementDiameter,
			double* u_dof,
			double* u_trial, 
			double* u_grad_trial, 
			double* u_test_dV, 
			double* u_grad_test_dV, 
			double* velocity,
			double* q_m,
			double* q_u,
			double* q_m_betaBDF,
			double* cfl,
			double* q_numDiff_u, 
			double* q_numDiff_u_last, 
			double* q_elementResidual_u, 
			int offset_u, int stride_u, 
			double* globalResidual,
			int nExteriorElementBoundaries_global,
			int* exteriorElementBoundariesArray,
			int* elementBoundaryElementsArray,
			int* elementBoundaryLocalElementBoundariesArray,
			double* u_trial_ext,
			double* u_grad_trial_ext,
			double* ebqe_velocity_ext,
			double* ebqe_n_ext,
			int* isDOFBoundary_u,
			double* ebqe_bc_u_ext,
			int* isFluxBoundary_u,
			double* ebqe_bc_flux_u_ext,
			double* u_test_dS_ext,
			double* ebqe_phi,double epsFact,
			double* ebqe_u,
			double* ebqe_flux)
{
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;

  //loop over elements to compute volume integrals and load them into element and global residual
  //
  //eN is the element index
  //eN_k is the quadrature point index for a scalar
  //eN_k_nSpace is the quadrature point index for a vector
  //eN_i is the element test function index
  //eN_j is the element trial function index
  //eN_k_j is the quadrature point index for a trial function
  //eN_k_i is the quadrature point index for a trial function
  double globalConservation=0.0;
  for(int eN=0;eN<nElements_global;eN++)
    {
      //declare local storage for element residual and initialize
      register double elementResidual_u[nDOF_test_element];
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_u[i]=0.0;
	}//i
      //loop over quadrature points and compute integrands
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  //compute indeces and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k,
	    eN_k_nSpace = eN_k*nSpace,
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double u=0.0,grad_u[nSpace],
	    m=0.0,dm=0.0,
	    f[nSpace],df[nSpace],
	    m_t=0.0,dm_t=0.0,
	    pdeResidual_u=0.0,
	    Lstar_u[nDOF_test_element],
	    subgridError_u=0.0,
	    tau=0.0,
            jac[nSpace*nSpace],
	    jacDet,
	    jacInv[nSpace*nSpace],
	    u_grad_trial[nDOF_trial_element*nSpace],
	    u_test_dV[nDOF_trial_element],
	    u_grad_test_dV[nDOF_test_element*nSpace],
	    dV,x,y,z,
	    G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv;
          // //
          // //compute solution and gradients at quadrature points
          // //
	  // u=0.0;
	  // for (int I=0;I<nSpace;I++)
	  //   {
	  //     grad_u[I]=0.0;
	  //   }
          // for (int j=0;j<nDOF_trial_element;j++)
          //   {
	  //     int eN_j=eN*nDOF_trial_element+j;
	  //     int eN_k_j=eN_k*nDOF_trial_element+j;
	  //     int eN_k_j_nSpace = eN_k_j*nSpace;
	  //     u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
	  //     for (int I=0;I<nSpace;I++)
	  //       {
	  //         grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
	  //       }
	  //   }
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
	  //get the trial function gradients
	  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	  //get the solution
	  ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	  //get the solution gradients
	  ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
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
          VOFV2_NAME::evaluateCoefficients(&velocity[eN_k_nSpace],
                                           u,
                                           m,
                                           dm,
                                           f,
                                           df);
          //
          //moving mesh
          //
          //omit for now
          //
          //calculate time derivative at quadrature points
          //
          ck.bdf(alphaBDF,
                 q_m_betaBDF[eN_k],
                 m,
                 dm,
                 m_t,
                 dm_t);
          //
          //calculate subgrid error (strong residual and adjoint)
          //
          //calculate strong residual
	  pdeResidual_u = ck.Mass_strong(m_t) + ck.Advection_strong(df,grad_u);
          //calculate adjoint
          for (int i=0;i<nDOF_test_element;i++)
            {
	      // register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
	      // Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);
	      register int i_nSpace = i*nSpace;
	      Lstar_u[i]  = ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
            }
          //calculate tau and tau*Res
          VOFV2_NAME::calculateSubgridError_tau(elementDiameter[eN],dm_t,df,cfl[eN_k],tau);
          subgridError_u = tau*pdeResidual_u;
          //
          //calcualte shock capturing diffusion
          //
          ck.calculateNumericalDiffusion(shockCapturingDiffusion,elementDiameter[eN],pdeResidual_u,grad_u,q_numDiff_u[eN_k]);
          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_k_i=eN_k*nDOF_test_element+i,
		eN_k_i_nSpace = eN_k_i*nSpace,
                i_nSpace=i*nSpace;

	      elementResidual_u[i] += ck.Mass_weak(m_t,u_test_dV[i]) + 
		ck.Advection_weak(f,&u_grad_test_dV[i_nSpace]) + 
 		ck.SubgridError(subgridError_u,Lstar_u[i]) + 
  		ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u,&u_grad_test_dV[i_nSpace]); 
	      
            }//i
	  //
	  //save solution for other models 
	  //
	  q_u[eN_k] = u;
	  q_m[eN_k] = m;
	}
      //
      //load element into global residual and save element residual
      //
      for(int i=0;i<nDOF_test_element;i++) 
        { 
          register int eN_i=eN*nDOF_test_element+i;
          
          q_elementResidual_u[eN_i] += elementResidual_u[i];
          globalConservation += elementResidual_u[i];
          globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
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
      const double eps=epsFact*elementDiameter[eN];
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
	    f_ext[nSpace],
	    df_ext[nSpace],
	    flux_ext=0.0,
	    bc_u_ext=0.0,
	    bc_grad_u_ext[nSpace],
	    bc_m_ext=0.0,
	    bc_dm_ext=0.0,
	    bc_f_ext[nSpace],
	    bc_df_ext[nSpace],
	    jac_ext[nSpace*nSpace],
	    jacDet_ext,
	    jacInv_ext[nSpace*nSpace],
	    boundaryJac[nSpace*(nSpace-1)],
	    metricTensor[(nSpace-1)*(nSpace-1)],
	    metricTensorDetSqrt,
	    dS,
            u_test_dS[nDOF_test_element],
	    u_grad_trial_trace[nDOF_trial_element*nSpace],
	    normal[3],x_ext,y_ext,z_ext,
	    G[nSpace*nSpace],G_dd_G,tr_G;
	  // 
	  //calculate the solution and gradients at quadrature points 
	  // 
	  // u_ext=0.0;
	  // for (int I=0;I<nSpace;I++)
	  //   {
	  //     grad_u_ext[I] = 0.0;
	  //     bc_grad_u_ext[I] = 0.0;/*mwf need way to evaluate this*/
	  //   }
	  // for (int j=0;j<nDOF_trial_element;j++) 
	  //   { 
	  //     int eN_j = eN*nDOF_trial_element+j;
	  //     int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
	  //     int ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	  //     u_ext += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial_ext[ebNE_kb_j]); 
	  //     for (int I=0;I<nSpace;I++)
	  //       {
	  //         grad_u_ext[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
	  //       } 
	  //   }
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
	  //calculate the pde coefficients using the solution and the boundary values for the solution 
	  // 
	  VOFV2_NAME::evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				 u_ext,
				 m_ext,
				 dm_ext,
				 f_ext,
				 df_ext);
          VOFV2_NAME::evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				 bc_u_ext,
				 bc_m_ext,
				 bc_dm_ext,
				 bc_f_ext,
				 bc_df_ext);    
	  //save for other models?
	  ebqe_u[ebNE_kb] = u_ext;
	  // 
	  //calculate the numerical fluxes 
	  // 
          VOFV2_NAME::exteriorNumericalAdvectiveFlux(isDOFBoundary_u[ebNE_kb],
                                                     isFluxBoundary_u[ebNE_kb],
                                                     normal,
                                                     bc_u_ext,
                                                     ebqe_bc_flux_u_ext[ebNE_kb],
                                                     u_ext,//smoothedHeaviside(eps,ebqe_phi[ebNE_kb]),
                                                     &ebqe_velocity_ext[ebNE_kb_nSpace],
                                                     flux_ext);
	  ebqe_flux[ebNE_kb] = flux_ext;
	  //
	  //update residuals
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;

	      elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS_ext[i]);
	    }//i
	}//kb
      //
      //update the element and global residual storage
      //
      for (int i=0;i<nDOF_test_element;i++)
	{
	  int eN_i = eN*nDOF_test_element+i;
	  q_elementResidual_u[eN_i] += elementResidual_u[i];
          globalConservation += elementResidual_u[i];
          globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
	}//i
    }//ebNE
  std::cout<<"VOFV2 global conservation============================================================="<<globalConservation<<std::endl;
}

//extern "C" void calculateJacobian_VOFV2(int nElements_global,
extern "C" void VOFV2_JAC (//element
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
			 double alphaBDF,
			 int lag_shockCapturing,/*mwf not used yet*/
			 double shockCapturingDiffusion,
			 int* u_l2g,
			 double* elementDiameter,
			 double* u_dof, 
			 double* u_trial, 
			 double* u_grad_trial, 
			 double* u_test_dV, 
			 double* u_grad_test_dV, 
			 double* velocity,
			 double* q_m_betaBDF, 
			 double* cfl,
			 double* q_numDiff_u_last, 
			 int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			 double* globalJacobian,
			 int nExteriorElementBoundaries_global,
			 int* exteriorElementBoundariesArray,
			 int* elementBoundaryElementsArray,
			 int* elementBoundaryLocalElementBoundariesArray,
			 double* u_trial_ext,
			 double* u_grad_trial_ext,
			 double* ebqe_velocity_ext,
			 double* ebqe_n_ext,
			 int* isDOFBoundary_u,
			 double* ebqe_bc_u_ext,
			 int* isFluxBoundary_u,
			 double* ebqe_bc_flux_u_ext,
			 double* u_test_dS_ext,
			 int* csrColumnOffsets_eb_u_u)
{
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;
  //
  //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
  //
  for(int eN=0;eN<nElements_global;eN++)
    {
      register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    elementJacobian_u_u[i][j]=0.0;
	  }
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
	    eN_k_nSpace = eN_k*nSpace,
	    eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	  //declare local storage
	  register double u=0.0,
	    grad_u[nSpace],
	    m=0.0,dm=0.0,
	    f[nSpace],df[nSpace],
	    m_t=0.0,dm_t=0.0,
	    dpdeResidual_u_u[nDOF_trial_element],
	    Lstar_u[nDOF_test_element],
	    dsubgridError_u_u[nDOF_trial_element],
	    tau=0.0,
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
	  // u=0.0;
	  // for (int I=0;I<nSpace;I++)
	  //   {
	  //     grad_u[I]=0.0;
	  //   }
          // for (int j=0;j<nDOF_trial_element;j++)
          //   {
	  //     int eN_j=eN*nDOF_trial_element+j;
	  //     int eN_k_j=eN_k*nDOF_trial_element+j;
	  //     int eN_k_j_nSpace = eN_k_j*nSpace;
              
          //     u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
	  //     for (int I=0;I<nSpace;I++)
	  //       {
	  //         grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
	  //       }
	  //   }
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
	  //get the trial function gradients
	  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	  //get the solution 	
	  ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	  //get the solution gradients
	  ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
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
          VOFV2_NAME::evaluateCoefficients(&velocity[eN_k_nSpace],
                                           u,
                                           m,
                                           dm,
                                           f,
                                           df);
          //
          //moving mesh
          //
          //omit for now
          //
          //calculate time derivatives
          //
          ck.bdf(alphaBDF,
                 q_m_betaBDF[eN_k],
                 m,
                 dm,
                 m_t,
                 dm_t);
          //
          //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
          //
          //calculate the adjoint times the test functions
          for (int i=0;i<nDOF_test_element;i++)
            {
	      // int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
	      // Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[eN_k_i_nSpace]);	      
	      register int i_nSpace = i*nSpace;
	      Lstar_u[i]=ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);	      
            }
          //calculate the Jacobian of strong residual
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
              int j_nSpace = j*nSpace;
	      dpdeResidual_u_u[j]= ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
		ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);
            }
          //tau and tau*Res
          VOFV2_NAME::calculateSubgridError_tau(elementDiameter[eN],
				      dm_t,
				      df,
				      cfl[eN_k],
				      tau);
          for(int j=0;j<nDOF_trial_element;j++)
            dsubgridError_u_u[j] = tau*dpdeResidual_u_u[j];
 	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_k_i=eN_k*nDOF_test_element+i;
	      int eN_k_i_nSpace=eN_k_i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int eN_k_j=eN_k*nDOF_trial_element+j;
		  int eN_k_j_nSpace = eN_k_j*nSpace;
		  int j_nSpace = j*nSpace;
		  int i_nSpace = i*nSpace;
		  elementJacobian_u_u[i][j] += ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) + 
		    ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]) +
		    ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) + 
		    ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]); 
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
	      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
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
  	    f_ext[nSpace],
  	    df_ext[nSpace],
  	    dflux_u_u_ext=0.0,
  	    bc_u_ext=0.0,
  	    bc_grad_u_ext[nSpace],
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
	    normal[3],x_ext,y_ext,z_ext,
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
	  dS = metricTensorDetSqrt*dS_ref[kb];
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
          VOFV2_NAME::evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
  				 u_ext,
  				 m_ext,
  				 dm_ext,
  				 f_ext,
  				 df_ext);
          VOFV2_NAME::evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
  				 bc_u_ext,
  				 bc_m_ext,
  				 bc_dm_ext,
  				 bc_f_ext,
  				 bc_df_ext);
  	  // 
  	  //calculate the numerical fluxes 
  	  // 
          VOFV2_NAME::exteriorNumericalAdvectiveFluxDerivative(isDOFBoundary_u[ebNE_kb],
                                                               isFluxBoundary_u[ebNE_kb],
                                                               normal,
                                                               &ebqe_velocity_ext[ebNE_kb_nSpace],
                                                               dflux_u_u_ext);
  	  //DoNothing for now
  	  //
  	  //calculate the flux jacobian
  	  //
  	  for (int j=0;j<nDOF_trial_element;j++)
  	    {
  	      register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
	      
  	      fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref[ebN_local_kb_j]);
  	    }//j
  	  //
  	  //update the global Jacobian from the flux Jacobian
  	  //
  	  for (int i=0;i<nDOF_test_element;i++)
  	    {
  	      register int eN_i = eN*nDOF_test_element+i,
  		ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
  	      for (int j=0;j<nDOF_trial_element;j++)
  		{
  		  register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;

  		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*u_test_dS_ext[ebNE_kb_i];
  		}//j
  	    }//i
  	}//kb
      }//ebNE
}//computeJacobian

