#include "MCorr2DQ.h"
#include <iostream>
#include <cassert>

//extern "C" void calculateResidual_MCorr(int nElements_global,
extern "C" void MCORR_RES(int nElements_global,
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
  using namespace MCORR_NAME;//MCorr;
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
	    eN_k_nSpace = eN_k*nSpace;
	  register double u=0.0,grad_u[nSpace],
	    r=0.0,dr=0.0;
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
          evaluateCoefficients_c(epsHeaviside,
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
		eN_k_i_nSpace = eN_k_i*nSpace;
	      
	      elementResidual_u[i] += Reaction_weak_c(r,u_test_dV[eN_k_i]) + 
		NumericalDiffusion_c(epsDiffusion*elementDiameter[eN],grad_u,&u_grad_test_dV[eN_k_i_nSpace]);
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

//extern "C" void calculateJacobian_MCorr(int nElements_global,
extern "C" void MCORR_JAC(int nElements_global,
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
  using namespace MCORR_NAME; //MCorr;
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
	    eN_k_nSpace = eN_k*nSpace; //index to a vector at a quadrature point

	  //declare local storage
	  register double u=0.0,
	    grad_u[nSpace],
	    r=0.0,dr=0.0;
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
          evaluateCoefficients_c(epsHeaviside,
				 epsDirac,
				 q_phi[eN_k],
				 q_H[eN_k],
				 u,
				 r,
				 dr);
 	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_k_i=eN_k*nDOF_test_element+i;
	      int eN_k_i_nSpace=eN_k_i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int eN_k_j=eN_k*nDOF_trial_element+j;
		  int eN_k_j_nSpace = eN_k_j*nSpace;
		  
		  elementJacobian_u_u[i][j] += ReactionJacobian_weak_c(dr,u_trial[eN_k_j],u_test_dV[eN_k_i]) + 
		    NumericalDiffusionJacobian_c(epsDiffusion*elementDiameter[eN],&u_grad_trial[eN_k_j_nSpace],&u_grad_test_dV[eN_k_i_nSpace]); 
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
