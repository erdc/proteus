#include <iostream>
#include "ElastoPlastic.h"
#include "CompKernel.h"
//plan
//1. Get it working with existing pre-compute approach
//a. Work out proper use of templates for inner dims (class for function + extern "C" calls)
//b. use linear elasticity first
//c. remember voigt rule is switched when translating old constitutive relations
//d. translate numerical stress flux functions
//2. Eliminate precomputations
//3. Potential improvements: indexing

extern "C" void calculateResidual_ElastoPlastic(int nElements_global,
						int* materialFlags,
						int nProperties_material,
						double* materialProperties,
						int* disp_l2g, 
						double* u_dof, 
						double* v_dof, 
						double* w_dof,
						double* disp_trial,
						double* disp_grad_trial,
						double* disp_test_dV,
						double* disp_grad_test_dV,
						double* force,
						int offset_u, 
						int offset_v, 
						int offset_w, 
						int stride_u, 
						int stride_v, 
						int stride_w, 
						double* globalResidual,
						int nExteriorElementBoundaries_global,
						int* exteriorElementBoundariesArray,
						int* elementBoundaryElementsArray,
						int* elementBoundaryLocalElementBoundariesArray,
						double* disp_trial_ext,
						double* disp_grad_trial_ext,
						double* disp_test_ext_dS,
						double* ebqe_force_ext,
						double* ebqe_normal_ext,
						int* isDispBoundary_u,
						int* isDispBoundary_v,
						int* isDispBoundary_w,
						int* isStressBoundary_u,
						int* isStressBoundary_v,
						int* isStressBoundary_w,
						double* ebqe_bc_u_ext,
						double* ebqe_bc_v_ext,
						double* ebqe_bc_w_ext,
						double* ebqe_bc_stress_u_ext,
						double* ebqe_bc_stress_v_ext,
						double* ebqe_bc_stress_w_ext,
						double* ebqe_penalty_ext)
{
  using namespace ElastoPlastic;
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
      register double elementResidual_u[nDOF_test_element],
	elementResidual_v[nDOF_test_element],
	elementResidual_w[nDOF_test_element];
      for (int i=0;i<nDOF_test_element;i++)
	{
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
	  register double u=0.0,v=0.0,w=0.0,grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],strain[nSym],stress[nSym],dstress[nSym2];
	  //
	  //mesh dependent info
	  //todo
	  //If we do this without pre-computing mesh-based info, then 
	  //Input: solution and mapping DOF, basis functions, gradients, Hessians for test, trial, and mapping spaces at quadrature points on reference element
	  //Compute physical points and Jacobian for space mapping
	  //Compute compute physical gradients and Hessians for test and trial functions
	  //
          //
          //compute solution and gradients at quadrature points
          //
	  u=0.0;v=0.0;w=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u[I]=0.0;grad_v[I]=0.0;grad_w[I]=0.0;
	    }
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_j=eN*nDOF_trial_element+j;
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
              u += valFromDOF_c(u_dof[disp_l2g[eN_j]],disp_trial[eN_k_j]);
              v += valFromDOF_c(v_dof[disp_l2g[eN_j]],disp_trial[eN_k_j]);
              w += valFromDOF_c(w_dof[disp_l2g[eN_j]],disp_trial[eN_k_j]);	      
	      for (int I=0;I<nSpace;I++)
		{
		  grad_u[I] += gradFromDOF_c(u_dof[disp_l2g[eN_j]],disp_grad_trial[eN_k_j_nSpace+I]);
		  grad_v[I] += gradFromDOF_c(v_dof[disp_l2g[eN_j]],disp_grad_trial[eN_k_j_nSpace+I]);
		  grad_w[I] += gradFromDOF_c(w_dof[disp_l2g[eN_j]],disp_grad_trial[eN_k_j_nSpace+I]);
		}
	    }
	  //Kinematic Voigt rule 
	  strain[0] = grad_u[0];
	  strain[1] = grad_v[1];
	  strain[2] = grad_w[2];
	  strain[3] = grad_v[2]+grad_w[1];
	  strain[4] = grad_u[2]+grad_w[0];
	  strain[5] = grad_u[1]+grad_v[0];
          //
          //calculate pde coefficients at quadrature points
          //
          evaluateCoefficients_c(materalProperties + nProperties_material*materialFlags[eN],//pointer to material array
				 strain,
				 stress,
				 dstress)
          //
          //moving mesh
          //
          //omit for now
          //
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_k_i=eN_k*nDOF_test_element+i,
		eN_k_i_nSpace = eN_k_i*nSpace;

	      elementResidual_u[i] += Stress_u_weak_c(stress,&disp_grad_test_dV[eN_k_i_nSpace]) + 
		Reaction_weak_c(force[eN_k],disp_test_dV[eN_k_i]);
	      elementResidual_v[i] += Stress_v_weak_c(stress,&disp_grad_test_dV[eN_k_i_nSpace]) + 
		Reaction_weak_c(force[eN_k],disp_test_dV[eN_k_i]);
	      elementResidual_w[i] += Stress_w_weak_c(stress,&disp_grad_test_dV[eN_k_i_nSpace]) + 
		Reaction_weak_c(force[eN_k],disp_test_dV[eN_k_i]);
            }//i
	}
      //
      //load element into global residual and save element residual
      //
      for(int i=0;i<nDOF_test_element;i++) 
        { 
          register int eN_i=eN*nDOF_test_element+i;
          
          q_elementResidual_u[eN_i]+=elementResidual_u[i];
          q_elementResidual_v[eN_i]+=elementResidual_v[i];
          q_elementResidual_w[eN_i]+=elementResidual_w[i];

          globalResidual[offset_u+stride_u*disp_l2g[eN_i]]+=elementResidual_u[i];
          globalResidual[offset_v+stride_v*disp_l2g[eN_i]]+=elementResidual_v[i];
          globalResidual[offset_w+stride_w*disp_l2g[eN_i]]+=elementResidual_w[i];
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
      register double elementResidual_u[nDOF_test_element],
	elementResidual_v[nDOF_test_element],
	elementResidual_w[nDOF_test_element];
      const double eps_rho = epsFact_rho*elementDiameter[eN],
      	eps_mu = epsFact_mu*elementDiameter[eN];
      for (int i=0;i<nDOF_test_element;i++)
	{
	  elementResidual_u[i]=0.0;
	  elementResidual_v[i]=0.0;
	  elementResidual_w[i]=0.0;
	}
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  register double u_ext=0.0,
	    v_ext=0.0,
	    w_ext=0.0,
	    grad_u_ext[nSpace],
	    grad_v_ext[nSpace],
	    grad_w_ext[nSpace],
	    bc_u_ext=0.0,
	    bc_v_ext=0.0,
	    bc_w_ext=0.0,
	    bc_stress_u_ext[nSpace],
	    bc_stress_v_ext[nSpace],
	    bc_stress_w_ext[nSpace];
	  // 
	  //calculate the solution and gradients at quadrature points 
	  // 
	 u_ext=0.0;v_ext=0.0;w_ext=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u_ext[I] = 0.0;
	      grad_v_ext[I] = 0.0;
	      grad_w_ext[I] = 0.0;
	      bc_grad_u_ext[I] = 0.0;
	      bc_grad_v_ext[I] = 0.0;
	      bc_grad_w_ext[I] = 0.0;
	    }
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      int eN_j = eN*nDOF_trial_element+j;
	      int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
	      int ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	      u_ext += valFromDOF_c(u_dof[disp_l2g[eN_j]],disp_trial_ext[ebNE_kb_j]); 
	      v_ext += valFromDOF_c(v_dof[disp_l2g[eN_j]],disp_trial_ext[ebNE_kb_j]); 
	      w_ext += valFromDOF_c(w_dof[disp_l2g[eN_j]],disp_trial_ext[ebNE_kb_j]); 
               
	      for (int I=0;I<nSpace;I++)
		{
		  grad_u_ext[I] += gradFromDOF_c(u_dof[disp_l2g[eN_j]],disp_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_v_ext[I] += gradFromDOF_c(v_dof[disp_l2g[eN_j]],disp_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_w_ext[I] += gradFromDOF_c(w_dof[disp_l2g[eN_j]],disp_grad_trial_ext[ebNE_kb_j_nSpace+I]);
		} 
	    }
	  //
	  //load the exterior trace (boundary condition) values
	  //
	  bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	  bc_v_ext = isDOFBoundary_v[ebNE_kb]*ebqe_bc_v_ext[ebNE_kb]+(1-isDOFBoundary_v[ebNE_kb])*v_ext;
	  bc_w_ext = isDOFBoundary_w[ebNE_kb]*ebqe_bc_w_ext[ebNE_kb]+(1-isDOFBoundary_w[ebNE_kb])*w_ext;
	  //Kinematic Voigt rule 
	  strain[0] = grad_u[0];
	  strain[1] = grad_v[1];
	  strain[2] = grad_w[2];
	  strain[3] = grad_v[2]+grad_w[1];
	  strain[4] = grad_u[2]+grad_w[0];
	  strain[5] = grad_u[1]+grad_v[0];
	  // 
	  //calculate the pde coefficients using the solution and the boundary values for the solution 
	  // 
          evaluateCoefficients_c(materalProperties + nProperties_material*materialFlags_ext[ebNE],//pointer to material array
				 strain,
				 stress,
				 dstress);
	  // 
	  //calculate the numerical fluxes 
	  // 
	  exteriorNumericalStressFlux(isDispBoundary_u[ebNE_kb],
				      isDispBoundary_v[ebNE_kb],
				      isDispBoundary_w[ebNE_kb],
				      isStressBoundary_u[ebNE_kb],
				      isStressBoundary_v[ebNE_kb],
				      isStressBoundary_w[ebNE_kb],
				      ebqe_n_ext[ebNE_kb_nSpace],
				      stress,
				      dstress,
				      bc_u_ext,
				      bc_v_ext,
				      bc_w_ext,
				      u_ext,
				      v_ext,
				      w_ext,
				      ebqe_penalty_ext[ebNE_kb],
				      ebqe_bc_stress_u_ext[ebNE_kb_nSpace],
				      ebqe_bc_stress_v_ext[ebNE_kb_nSpace],
				      ebqe_bc_stress_w_ext[ebNE_kb_nSpace],
				      stressFlux_u,
				      stressFlux_v,
				      stressFlux_w)

	  //
	  //update residuals
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
	      elementResidual_u[i] += ExteriorElementBoundaryFlux_c(stressFlux_u,disp_test_dS_ext[ebNE_kb_i]);
	      elementResidual_v[i] += ExteriorElementBoundaryFlux_c(stressFlux_v,disp_test_dS_ext[ebNE_kb_i]);	       
	      elementResidual_w[i] += ExteriorElementBoundaryFlux_c(stressFlux_w,disp_test_dS_ext[ebNE_kb_i]);
	    }//i
	}//kb
      //
      //update the element and global residual storage
      //
      for (int i=0;i<nDOF_test_element;i++)
	{
	  int eN_i = eN*nDOF_test_element+i;
	  globalResidual[offset_u+stride_u*disp_l2g[eN_i]]+=elementResidual_u[i];
	  globalResidual[offset_v+stride_v*disp_l2g[eN_i]]+=elementResidual_v[i];
	  globalResidual[offset_w+stride_w*disp_l2g[eN_i]]+=elementResidual_w[i];
	}//i
    }//ebNE
}

extern "C" void calculateJacobian_ElastoPlastic(int nElements_global,
						int nProperties_material,
						double* materalProperties,
						int* disp_l2g,
						double* u_dof, 
						double* v_dof, 
						double* w_dof,
						double* disp_trial,
						double* disp_grad_trial,
						double* disp_test_dV, 
						double* disp_grad_test_dV,
						double* force,
						double* normal,
						int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
						int* csrRowIndeces_u_v,int* csrColumnOffsets_u_v,
						int* csrRowIndeces_u_w,int* csrColumnOffsets_u_w,
						int* csrRowIndeces_v_u,int* csrColumnOffsets_v_u,
						int* csrRowIndeces_v_v,int* csrColumnOffsets_v_v,
						int* csrRowIndeces_v_w,int* csrColumnOffsets_v_w,
						int* csrRowIndeces_w_u,int* csrColumnOffsets_w_u,
						int* csrRowIndeces_w_v,int* csrColumnOffsets_w_v,
						int* csrRowIndeces_w_w,int* csrColumnOffsets_w_w,
						double* globalJacobian,
						int nExteriorElementBoundaries_global,
						int* exteriorElementBoundariesArray,
						int* elementBoundaryElementsArray,
						int* elementBoundaryLocalElementBoundariesArray,
						double* disp_trial_ext,
						double* disp_grad_trial_ext,
						double* ebqe_force_ext,
						double* ebqe_normal_ext,
						int* isDispBoundary_u,
						int* isDispBoundary_v,
						int* isDispBoundary_w,
						int* isStressBoundary_u,
						int* isStressBoundary_v,
						int* isStressBoundary_w,
						double* ebqe_bc_stress_u_ext,
						double* ebqe_bc_stress_v_ext,
						double* ebqe_bc_stress_w_ext,
						double* ebqe_penalty_ext,
						double* disp_test_dS_ext)
{
  using namespace ElastoPlastic;
  //
  //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
  //
  for(int eN=0;eN<nElements_global;eN++)
    {
      register double elementJacobian_u_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_u_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_v_w[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_u[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_v[nDOF_test_element][nDOF_trial_element],
	elementJacobian_w_w[nDOF_test_element][nDOF_trial_element];
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    elementJacobian_u_u[i][j]=0.0;
	    elementJacobian_u_v[i][j]=0.0;
	    elementJacobian_u_w[i][j]=0.0;
	    elementJacobian_v_u[i][j]=0.0;
	    elementJacobian_v_v[i][j]=0.0;
	    elementJacobian_v_w[i][j]=0.0;
	    elementJacobian_w_u[i][j]=0.0;
	    elementJacobian_w_v[i][j]=0.0;
	    elementJacobian_w_w[i][j]=0.0;
	  }
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
	    eN_k_nSpace = eN_k*nSpace; //index to a vector at a quadrature point

	  //declare local storage
	  register double u=0.0,v=0.0,w=0.0,
	    grad_u[nSpace],grad_v[nSpace],grad_w[nSpace];
          //
          //calculate solution and gradients at quadrature points
          //
	  u=0.0;v=0.0;w=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u[I]=0.0;grad_v[I]=0.0;grad_w[I]=0.0;
	    }
          for (int j=0;j<nDOF_trial_element;j++)
            {
	      int eN_j=eN*nDOF_trial_element+j;
	      int eN_k_j=eN_k*nDOF_trial_element+j;
	      int eN_k_j_nSpace = eN_k_j*nSpace;
              u += valFromDOF_c(u_dof[disp_l2g[eN_j]],disp_trial[eN_k_j]);
              v += valFromDOF_c(v_dof[disp_l2g[eN_j]],disp_trial[eN_k_j]);
              w += valFromDOF_c(w_dof[disp_l2g[eN_j]],disp_trial[eN_k_j]);
	      for (int I=0;I<nSpace;I++)
		{
		  grad_u[I] += gradFromDOF_c(u_dof[disp_l2g[eN_j]],disp_grad_trial[eN_k_j_nSpace+I]);
		  grad_v[I] += gradFromDOF_c(v_dof[disp_l2g[eN_j]],disp_grad_trial[eN_k_j_nSpace+I]);
		  grad_w[I] += gradFromDOF_c(w_dof[disp_l2g[eN_j]],disp_grad_trial[eN_k_j_nSpace+I]);
		}
	    }
          //
          //calculate pde coefficients and derivatives at quadrature points
          //
          evaluateCoefficients_c(materalProperties + nProperties_material*materialFlags[eN],//pointer to material array
				 strain,
				 stress,
				 dstress);
	  //convert dstress to dStress_u_u, etc.
	  //
          //
          //moving mesh
          //
          //omit for now
          //
  	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_k_i=eN_k*nDOF_test_element+i;
	      int eN_k_i_nSpace=eN_k_i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int eN_k_j=eN_k*nDOF_trial_element+j;
		  int eN_k_j_nSpace = eN_k_j*nSpace;
		  elementJacobian_u_u[i][j] += stressJacobian_weak_c(dStress_u_u,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_u_v[i][j] += stressJacobian_weak_c(dStress_u_v,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_u_w[i][j] += stressJacobian_weak_c(dStress_u_w,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_v_u[i][j] += stressJacobian_weak_c(dStress_v_u,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_v_v[i][j] += stressJacobian_weak_c(dStress_v_v,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_v_w[i][j] += stressJacobian_weak_c(dStress_v_w,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_w_u[i][j] += stressJacobian_weak_c(dStress_w_u,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_w_v[i][j] += stressJacobian_weak_c(dStress_w_v,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
		  elementJacobian_w_w[i][j] += stressJacobian_weak_c(dStress_w_w,&disp_grad_trial[eN_k_j_nSpace],&disp_grad_test_dV[eN_k_i_nSpace]);
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
	      globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_u_v[eN_i_j]] += elementJacobian_u_v[i][j];
	      globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_u_w[eN_i_j]] += elementJacobian_u_w[i][j];

	      globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_v_u[eN_i_j]] += elementJacobian_v_u[i][j];
	      globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_v_v[eN_i_j]] += elementJacobian_v_v[i][j];
	      globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_v_w[eN_i_j]] += elementJacobian_v_w[i][j];

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
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	{ 
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;

	  register double u_ext=0.0,
	    v_ext=0.0,
	    w_ext=0.0,
	    grad_u_ext[nSpace],
	    grad_v_ext[nSpace],
	    grad_w_ext[nSpace],
	    bc_u_ext=0.0,
	    bc_v_ext=0.0,
	    bc_w_ext=0.0;
	  //
	  //calculate mesh dependent FEM info
	  //
	  //todo
	  // 
	  //calculate the solution and gradients at quadrature points 
	  // 
	  u_ext=0.0;v_ext=0.0;w_ext=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u_ext[I] = 0.0;
	      grad_v_ext[I] = 0.0;
	      grad_w_ext[I] = 0.0;
	      bc_grad_u_ext[I] = 0.0;
	      bc_grad_v_ext[I] = 0.0;
	      bc_grad_w_ext[I] = 0.0;
	    }
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      register int eN_j = eN*nDOF_trial_element+j,
		ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
		ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
	      u_ext += valFromDOF_c(u_dof[disp_l2g[eN_j]],disp_trial_ext[ebNE_kb_j]); 
	      v_ext += valFromDOF_c(v_dof[disp_l2g[eN_j]],disp_trial_ext[ebNE_kb_j]); 
	      w_ext += valFromDOF_c(w_dof[disp_l2g[eN_j]],disp_trial_ext[ebNE_kb_j]); 
               
	      for (int I=0;I<nSpace;I++)
		{
		  grad_u_ext[I] += gradFromDOF_c(u_dof[disp_l2g[eN_j]],disp_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_v_ext[I] += gradFromDOF_c(v_dof[disp_l2g[eN_j]],disp_grad_trial_ext[ebNE_kb_j_nSpace+I]); 
		  grad_w_ext[I] += gradFromDOF_c(w_dof[disp_l2g[eN_j]],disp_grad_trial_ext[ebNE_kb_j_nSpace+I]);
		} 
	    }
	  //
	  //load the boundary values
	  //
	  bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	  bc_v_ext = isDOFBoundary_v[ebNE_kb]*ebqe_bc_v_ext[ebNE_kb]+(1-isDOFBoundary_v[ebNE_kb])*v_ext;
	  bc_w_ext = isDOFBoundary_w[ebNE_kb]*ebqe_bc_w_ext[ebNE_kb]+(1-isDOFBoundary_w[ebNE_kb])*w_ext;
	  // 
	  //calculate the internal and external trace of the pde coefficients 
	  // 
	  evaluateCoefficinets_c(...);
	  evaluateCoefficients_c(...);
	  // 
	  //calculate the numerical fluxes 
	  // 
	  exteriorNumericalStressFluxJacobian(isDispBoundary_u[ebNE_kb],
					      isDispBoundary_v[ebNE_kb],
					      isDispBoundary_w[ebNE_kb],
					      isStressBoundary_u[ebNE_kb],
					      isStressBoundary_v[ebNE_kb],
					      isStressBoundary_w[ebNE_kb],
					      ebqe_n_ext[ebNE_kb_nSpace],
					      stress,
					      dstress,
					      bc_u_ext,
					      bc_v_ext,
					      bc_w_ext,
					      u_ext,
					      v_ext,
					      w_ext,
					      ebqe_penalty_ext[ebNE_kb],
					      ebqe_bc_stress_u_ext[ebNE_kb_nSpace],
					      ebqe_bc_stress_v_ext[ebNE_kb_nSpace],
					      ebqe_bc_stress_w_ext[ebNE_kb_nSpace],
					      dstressFlux_u_u,
					      dstressFlux_u_v,
					      dstressFlux_u_w,
					      dstressFlux_v_u,
					      dstressFlux_v_v,
					      dstressFlux_v_w,
					      dstressFlux_w_u,
					      dstressFlux_w_v,
					      dstressFlux_w_w);
	  //
	  //calculate the flux jacobian
	  //
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
		ebNE_kb_j_nSpace = ebNE_kb_j*nSpace;

	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_u_u,disp_trial_ext[ebNE_kb_j]);
	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_u_v,disp_trial_ext[ebNE_kb_j]);
	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_u_w,disp_trial_ext[ebNE_kb_j]);

	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_v_u,disp_trial_ext[ebNE_kb_j]);
	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_v_v,disp_trial_ext[ebNE_kb_j]);
	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_v_w,disp_trial_ext[ebNE_kb_j]);

	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_w_u,disp_trial_ext[ebNE_kb_j]);
	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_w_v,disp_trial_ext[ebNE_kb_j]);
	      fluxJacobian_u_u[j]=ExteriorNumericalStressFluxJacobian_c(dstressFlux_w_w,disp_trial_ext[ebNE_kb_j]);
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
		  
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*disp_test_dS_ext[ebNE_kb_i];
		  globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_eb_u_v[ebN_i_j]] += fluxJacobian_u_v[j]*disp_test_dS_ext[ebNE_kb_i];
		  globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_eb_u_w[ebN_i_j]] += fluxJacobian_u_w[j]*disp_test_dS_ext[ebNE_kb_i];
		   
		  globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_eb_v_u[ebN_i_j]] += fluxJacobian_v_u[j]*disp_test_dS_ext[ebNE_kb_i];
		  globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_eb_v_v[ebN_i_j]] += fluxJacobian_v_v[j]*disp_test_dS_ext[ebNE_kb_i];
		  globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_eb_v_w[ebN_i_j]] += fluxJacobian_v_w[j]*disp_test_dS_ext[ebNE_kb_i];
		   
		  globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_eb_w_u[ebN_i_j]] += fluxJacobian_w_u[j]*disp_test_dS_ext[ebNE_kb_i];
		  globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_eb_w_v[ebN_i_j]] += fluxJacobian_w_v[j]*disp_test_dS_ext[ebNE_kb_i];
		  globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_eb_w_w[ebN_i_j]] += fluxJacobian_w_w[j]*disp_test_dS_ext[ebNE_kb_i];
		}//j
	    }//i
	}//kb
    }//ebNE
}//computeJacobian
