#include "MCorrV2.h"
#include <iostream>

//extern "C" void calculateResidual_MCorrV2(int nElements_global,
extern "C" void MCORRV2_RES(//element
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
                            double epsFactHeaviside,
                            double epsFactDirac,
                            double epsFactDiffusion,
                            int* u_l2g, 
                            double* elementDiameter,
                            double* u_dof,
                            double* u_trial, 
                            double* u_grad_trial, 
                            double* u_test_dV, 
                            double* u_grad_test_dV, 
                            double* q_phi,
                            double* q_H,
                            double* q_u,
                            double* q_r,
                            int offset_u, int stride_u, 
                            double* globalResidual)
{
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;
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
      register double elementResidual_u[nDOF_test_element];
      const double epsHeaviside=epsFactHeaviside*elementDiameter[eN],
	epsDirac=epsFactDirac*elementDiameter[eN],
	epsDiffusion=epsFactDiffusion*elementDiameter[eN];
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
	    r=0.0,dr=0.0,
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
          MCORRV2_NAME::evaluateCoefficients(epsHeaviside,
                                             epsDirac,
                                             q_phi[eN_k],
                                             q_H[eN_k],
                                             u,
                                             r,
                                             dr);
          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_k_i=eN_k*nDOF_test_element+i,
		eN_k_i_nSpace = eN_k_i*nSpace,
                i_nSpace=i*nSpace;
	      
	      elementResidual_u[i] += ck.Reaction_weak(r,u_test_dV[i]) + 
		ck.NumericalDiffusion(epsDiffusion*elementDiameter[eN],grad_u,&u_grad_test_dV[i_nSpace]);
            }//i
	  //
	  //save momentum for time history and velocity for subgrid error
	  //save solution for other models 
	  //
	  q_r[eN_k] = r;
	  q_u[eN_k] = u;
	}
      //
      //load element into global residual and save element residual
      //
      for(int i=0;i<nDOF_test_element;i++) 
        { 
          register int eN_i=eN*nDOF_test_element+i;
          
          globalResidual[offset_u+stride_u*u_l2g[eN_i]]+=elementResidual_u[i];
        }//i
    }//elements
}

//extern "C" void calculateJacobian_MCorrV2(int nElements_global,
extern "C" void MCORRV2_JAC(//element
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
                          double epsFactHeaviside,
                          double epsFactDirac,
                          double epsFactDiffusion,
                          int* u_l2g,
                          double* elementDiameter,
                          double* u_dof, 
                          double* u_trial, 
                          double* u_grad_trial, 
                          double* u_test_dV, 
                          double* u_grad_test_dV, 
                          double* q_phi,
                          double* q_H,
                          int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
                          double* globalJacobian)
{
  CompKernel<nSpace,nDOF_mesh_trial_element,nDOF_trial_element,nDOF_test_element> ck;
  //
  //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
  //
  for(int eN=0;eN<nElements_global;eN++)
    {
      register double  elementJacobian_u_u[nDOF_test_element][nDOF_trial_element];
      const double epsHeaviside=epsFactHeaviside*elementDiameter[eN],
	epsDirac=epsFactDirac*elementDiameter[eN],
	epsDiffusion=epsFactDiffusion*elementDiameter[eN];
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
	    r=0.0,dr=0.0,
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
          MCORRV2_NAME::evaluateCoefficients(epsHeaviside,
                                               epsDirac,
                                               q_phi[eN_k],
                                               q_H[eN_k],
                                               u,
                                               r,
                                               dr);
 	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_k_i=eN_k*nDOF_test_element+i;
	      int eN_k_i_nSpace=eN_k_i*nSpace,
                i_nSpace=i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int eN_k_j=eN_k*nDOF_trial_element+j;
		  int eN_k_j_nSpace = eN_k_j*nSpace,
                    j_nSpace = j*nSpace;
		  
		  elementJacobian_u_u[i][j] += ck.ReactionJacobian_weak(dr,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) + 
		    ck.NumericalDiffusionJacobian(epsDiffusion*elementDiameter[eN],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]); 
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
}//computeJacobian
