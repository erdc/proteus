#include "VOF.h"
#include <iostream>

//extern "C" void calculateResidual_VOF(int nElements_global,
extern "C" void VOF_RES(int nElements_global,
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
			double* q_cfl,
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
  using namespace VOF_NAME;
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
	    eN_k_nSpace = eN_k*nSpace;
	  register double u=0.0,grad_u[nSpace],
	    m=0.0,dm=0.0,
	    f[nSpace],df[nSpace],
	    m_t=0.0,dm_t=0.0,
	    pdeResidual_u=0.0,
	    Lstar_u[nDOF_test_element],
	    subgridError_u=0.0,
	    tau=0.0;
          //
          //compute solution and gradients at quadrature points
          //
	  u=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u[I]=0.0;
	    }
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_j=eN*nDOF_trial_element+j;
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
	      u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
	      for (int I=0;I<nSpace;I++)
		{
		  grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
		}
	    }
          //
          //calculate pde coefficients at quadrature points
          //
          evaluateCoefficients_c(&velocity[eN_k_nSpace],
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
          bdf_c(alphaBDF,
		q_m_betaBDF[eN_k],
		m,
		dm,
		m_t,
		dm_t);
          //
          //calculate subgrid error (strong residual and adjoint)
          //
          //calculate strong residual
	  pdeResidual_u = Mass_strong_c(m_t) + Advection_strong_c(df,grad_u);
          //calculate adjoint
          for (int i=0;i<nDOF_test_element;i++)
            {
	      register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
	      Lstar_u[i]  = Advection_adjoint_c(df,&u_grad_test_dV[eN_k_i_nSpace]);
            }
          //calculate tau and tau*Res
          calculateSubgridError_tau_c(elementDiameter[eN],dm_t,df,q_cfl[eN_k],tau);
          calculateSubgridError_tauRes_c(tau,pdeResidual_u,subgridError_u);
          //
          //calcualte shock capturing diffusion
          //
          calculateNumericalDiffusion_c(shockCapturingDiffusion,elementDiameter[eN],pdeResidual_u,grad_u,q_numDiff_u[eN_k]);
          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_k_i=eN_k*nDOF_test_element+i,
		eN_k_i_nSpace = eN_k_i*nSpace;

	      elementResidual_u[i] += Mass_weak_c(m_t,u_test_dV[eN_k_i]) + 
		Advection_weak_c(f,&u_grad_test_dV[eN_k_i_nSpace]) + 
 		SubgridError_c(subgridError_u,Lstar_u[i]) + 
  		NumericalDiffusion_c(q_numDiff_u_last[eN_k],grad_u,&u_grad_test_dV[eN_k_i_nSpace]); 
	      
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
	eN  = elementBoundaryElementsArray[ebN*2+0];
      register double elementResidual_u[nDOF_test_element];
      const double eps=epsFact*elementDiameter[eN];
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_u[i]=0.0;
	}
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
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
	    bc_df_ext[nSpace];
	  // 
	  //calculate the solution and gradients at quadrature points 
	  // 
	  u_ext=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u_ext[I] = 0.0;
	      bc_grad_u_ext[I] = 0.0;/*mwf need way to evaluate this*/
	    }
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      int eN_j = eN*nDOF_trial_element+j;
	      int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
	      int ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	      u_ext += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial_ext[ebNE_kb_j]); 
	      for (int I=0;I<nSpace;I++)
		{
		  grad_u_ext[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		} 
	    }
	  //
	  //load the boundary values
	  //
	  bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;

	  // 
	  //calculate the pde coefficients using the solution and the boundary values for the solution 
	  // 
	  evaluateCoefficients_c(&ebqe_velocity_ext[ebNE_kb_nSpace],
				 u_ext,
				 m_ext,
				 dm_ext,
				 f_ext,
				 df_ext);
	  evaluateCoefficients_c(&ebqe_velocity_ext[ebNE_kb_nSpace],
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
	  exteriorNumericalAdvectiveFlux_c(isDOFBoundary_u[ebNE_kb],
					   isFluxBoundary_u[ebNE_kb],
					   &ebqe_n_ext[ebNE_kb_nSpace],
					   bc_u_ext,
					   ebqe_bc_flux_u_ext[ebNE_kb],
					   smoothedHeaviside(eps,ebqe_phi[ebNE_kb]),
					   &ebqe_velocity_ext[ebNE_kb_nSpace],
					   flux_ext);
	  ebqe_flux[ebNE_kb] = flux_ext;
	  //
	  //update residuals
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;

	      elementResidual_u[i] += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
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
  std::cout<<"VOF global conservation============================================================="<<globalConservation<<std::endl;
}

//extern "C" void calculateJacobian_VOF(int nElements_global,
extern "C" void VOF_JAC (int nElements_global,
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
			 double* q_cfl,
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
  using namespace VOF_NAME;
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
	    eN_k_nSpace = eN_k*nSpace; //index to a vector at a quadrature point

	  //declare local storage
	  register double u=0.0,
	    grad_u[nSpace],
	    m=0.0,dm=0.0,
	    f[nSpace],df[nSpace],
	    m_t=0.0,dm_t=0.0,
	    dpdeResidual_u_u[nDOF_trial_element],
	    Lstar_u[nDOF_test_element],
	    dsubgridError_u_u[nDOF_trial_element],
	    tau=0.0;
          //
          //calculate solution and gradients at quadrature points
          //
	  u=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u[I]=0.0;
	    }
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_j=eN*nDOF_trial_element+j;
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
              
              u += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial[eN_k_j]);
	      for (int I=0;I<nSpace;I++)
		{
		  grad_u[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial[eN_k_j_nSpace+I]);
		}
	    }
          //
          //calculate pde coefficients and derivatives at quadrature points
          //
          evaluateCoefficients_c(&velocity[eN_k_nSpace],
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
          bdf_c(alphaBDF,
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
	      int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;

	      Lstar_u[i]=Advection_adjoint_c(df,&u_grad_test_dV[eN_k_i_nSpace]);	      
            }
          //calculate the Jacobian of strong residual
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
	      dpdeResidual_u_u[j]= MassJacobian_strong_c(dm_t,u_trial[eN_k_j]) +
		AdvectionJacobian_strong_c(df,&u_grad_trial[eN_k_j_nSpace]);
            }
          //tau and tau*Res
          calculateSubgridError_tau_c(elementDiameter[eN],
				      dm_t,
				      df,
				      q_cfl[eN_k],
				      tau);
          calculateSubgridErrorDerivatives_tauRes_c(tau,
						    dpdeResidual_u_u,
						    dsubgridError_u_u);
 	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_k_i=eN_k*nDOF_test_element+i;
	      int eN_k_i_nSpace=eN_k_i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int eN_k_j=eN_k*nDOF_trial_element+j;
		  int eN_k_j_nSpace = eN_k_j*nSpace;
		  
		  elementJacobian_u_u[i][j] += MassJacobian_weak_c(dm_t,u_trial[eN_k_j],u_test_dV[eN_k_i]) + 
		    AdvectionJacobian_weak_c(df,u_trial[eN_k_j],&u_grad_test_dV[eN_k_i_nSpace]) +
		    SubgridErrorJacobian_c(dsubgridError_u_u[j],Lstar_u[i]) + 
		    NumericalDiffusionJacobian_c(q_numDiff_u_last[eN_k],&u_grad_trial[eN_k_j_nSpace],&u_grad_test_dV[eN_k_i_nSpace]); 		  
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
  //   //
  //   //loop over exterior element boundaries to compute the surface integrals and load them into the global Jacobian
  //   //
  //   for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
  //     { 
  //       register int ebN = exteriorElementBoundariesArray[ebNE]; 
  //       register int eN  = elementBoundaryElementsArray[ebN*2+0];
  //       for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
  // 	{ 
  // 	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
  // 	    ebNE_kb_nSpace = ebNE_kb*nSpace;

  // 	  register double u_ext=0.0,
  // 	    grad_u_ext[nSpace],
  // 	    m_ext=0.0,
  // 	    dm_ext=0.0,
  // 	    f_ext[nSpace],
  // 	    df_ext[nSpace],
  // 	    dflux_u_u_ext=0.0,
  // 	    bc_u_ext=0.0,
  // 	    bc_grad_u_ext[nSpace],
  // 	    bc_m_ext=0.0,
  // 	    bc_dm_ext=0.0,
  // 	    bc_f_ext[nSpace],
  // 	    bc_df_ext[nSpace],
  // 	    fluxJacobian_u_u[nDOF_trial_element];
  // 	  // 
  // 	  //calculate the solution and gradients at quadrature points 
  // 	  // 
  // 	  u_ext=0.0;
  // 	  for (int I=0;I<nSpace;I++)
  // 	    {
  // 	      grad_u_ext[I] = 0.0;
  // 	      bc_grad_u_ext[I] = 0.0;
  // 	    }
  // 	  for (int j=0;j<nDOF_trial_element;j++) 
  // 	    { 
  // 	      register int eN_j = eN*nDOF_trial_element+j,
  // 		ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
  // 		ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
  // 	      u_ext += valFromDOF_c(u_dof[u_l2g[eN_j]],u_trial_ext[ebNE_kb_j]); 
	                     
  // 	      for (int I=0;I<nSpace;I++)
  // 		{
  // 		  grad_u_ext[I] += gradFromDOF_c(u_dof[u_l2g[eN_j]],u_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
  // 		} 
  // 	    }
  // 	  //
  // 	  //load the boundary values
  // 	  //
  // 	  bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
  // 	  // 
  // 	  //calculate the internal and external trace of the pde coefficients 
  // 	  // 
  // 	  evaluateCoefficients_c(&ebqe_velocity_ext[ebNE_kb_nSpace],
  // 				 u_ext,
  // 				 m_ext,
  // 				 dm_ext,
  // 				 f_ext,
  // 				 df_ext);
  // 	  evaluateCoefficients_c(&ebqe_velocity_ext[ebNE_kb_nSpace],
  // 				 bc_u_ext,
  // 				 bc_m_ext,
  // 				 bc_dm_ext,
  // 				 bc_f_ext,
  // 				 bc_df_ext);
  // 	  // 
  // 	  //calculate the numerical fluxes 
  // 	  // 
  // 	  exteriorNumericalAdvectiveFluxDerivative_c(isDOFBoundary_u[ebNE_kb],
  // 						     isFluxBoundary_u[ebNE_kb],
  // 						     &ebqe_n_ext[ebNE_kb_nSpace],
  // 						     &ebqe_velocity_ext[ebNE_kb_nSpace],
  // 						     dflux_u_u_ext);
  // 	  //DoNothing for now
  // 	  //
  // 	  //calculate the flux jacobian
  // 	  //
  // 	  for (int j=0;j<nDOF_trial_element;j++)
  // 	    {
  // 	      register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
	      
  // 	      fluxJacobian_u_u[j]=ExteriorNumericalAdvectiveFluxJacobian_c(dflux_u_u_ext,u_trial_ext[ebNE_kb_j]);
  // 	    }//j
  // 	  //
  // 	  //update the global Jacobian from the flux Jacobian
  // 	  //
  // 	  for (int i=0;i<nDOF_test_element;i++)
  // 	    {
  // 	      register int eN_i = eN*nDOF_test_element+i,
  // 		ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
  // 	      for (int j=0;j<nDOF_trial_element;j++)
  // 		{
  // 		  register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;

  // 		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*u_test_dS_ext[ebNE_kb_i];
  // 		}//j
  // 	    }//i
  // 	}//kb
  //     }//ebNE
}//computeJacobian

