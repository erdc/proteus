#ifndef PRES_H
#define PRES_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

namespace proteus
{
  class cppPres_base
  {
  public:
    virtual ~cppPres_base(){}
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
				   double* q_u,
				   double* q_grad_u,
				   double* q_p_last,
				   double* q_p_inc,
				   double* q_massFlux,
				   double* ebqe_massFlux,
				   double* ebqe_u,
				   double* ebqe_grad_u,
				   int offset_u, int stride_u, 
				   double* globalResidual,			   
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray)=0;
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
				   int nElements_global,
				   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				   double* globalJacobian)=0;
  };
  
  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class cppPres : public cppPres_base
  {
  public:
    CompKernelType ck;
    cppPres():ck()
    {}
    inline
      void evaluateCoefficients(const double& p,
				const double& p_last,
				const double& p_inc,
				const double massFlux[nSpace],
				double f[nSpace],
				double& r)
    {
      for (int I=0;I<nSpace;I++)
        f[I] = massFlux[I];
      r = p - p_last - p_inc;
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
					 double* q_u,
					 double* q_grad_u,
					 double* q_p_last,
					 double* q_p_inc,
					 double* q_massFlux,
					 double* ebqe_massFlux,
					 double* ebqe_u,
					 double* ebqe_grad_u,
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
      double epsHeaviside,epsDirac,epsDiffusion,norm;
      //loop over quadrature points and compute integrands
      for  (int k=0;k<nQuadraturePoints_element;k++)
	{
	  //compute indeces and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k,
	    eN_k_nSpace = eN_k*nSpace;
	    //eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double u=0.0,grad_u[nSpace],
	    f[nSpace],
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
	  ck.valFromElementDOF(element_u,&u_trial_ref[k*nDOF_trial_element],u);
	  //get the solution gradients
	  ck.gradFromElementDOF(element_u,u_grad_trial,grad_u);
	  //precalculate test function products with integration weights
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
	    }
	  //
	  //calculate pde coefficients at quadrature points
	  //
	  evaluateCoefficients(u,//u is p
			       q_p_last[eN_k],
			       q_p_inc[eN_k],
			       &q_massFlux[eN_k_nSpace],
			       f,
			       r);
	  // 
	  //update element residual 
	  // 
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int  i_nSpace=i*nSpace;	      
	      elementResidual_u[i] += ck.Advection_weak(f,&u_grad_test_dV[i_nSpace]) +
		ck.Reaction_weak(r,u_test_dV[i]);
	    }//i
	  q_u[eN_k] = u;
	  for(int I=0;I<nSpace;I++)
	    q_grad_u[eN_k_nSpace+I] = grad_u[I];
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
			   double* q_u,
			   double* q_grad_u,
			   double* q_p_last,
			   double* q_p_inc,
			   double* q_massFlux,
			   double* ebqe_massFlux,
			   double* ebqe_u,
			   double* ebqe_grad_u,
			   int offset_u, int stride_u, 
			   double* globalResidual,			   
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray)
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
				   q_u,
				   q_grad_u,
				   q_p_last,
				   q_p_inc,
				   q_massFlux,
				   ebqe_massFlux,
				   ebqe_u,
				   ebqe_grad_u,
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
      //loop over exterior element boundaries to calculate levelset gradient
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE], 
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
	    //eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
	  double element_u[nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i=eN*nDOF_test_element+i;
	      element_u[i] = u_dof[u_l2g[eN_i]];
              elementResidual_u[i] = 0.0;
	    }//i
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	      register double penalty=0.0,
                u_ext=0.0,
                bc_u_ext=0.0,
                adv_flux_ext=0.0,
                diff_flux_ext=0.0,
                a_ext,
                f_ext[nSpace],
		grad_u_ext[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		u_grad_trial_trace[nDOF_trial_element*nSpace],
		u_grad_test_dS[nDOF_test_element*nSpace],
		normal[nSpace],x_ext,y_ext,z_ext,
		G[nSpace*nSpace],G_dd_G,tr_G;
	      // 
	      //calculate the solution and gradients at quadrature points 
	      // 
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
	      ck.calculateGScale(G,normal,penalty);
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&u_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,u_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromElementDOF(element_u,&u_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.gradFromElementDOF(element_u,u_grad_trial_trace,grad_u_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dS[j] = u_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		}
	      //
	      //load the boundary values
	      //
	      ebqe_u[ebNE_kb] = u_ext;
	      for (int I=0;I<nSpace;I++)
		ebqe_grad_u[ebNE_kb_nSpace+I] = grad_u_ext[I];
	      adv_flux_ext=0.0;
	      for (int I=0;I<nSpace;I++)
		adv_flux_ext += ebqe_massFlux[ebNE_kb_nSpace+I]*normal[I];
	      //
	      //update residuals
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(adv_flux_ext,u_test_dS[i]);
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
					 double* elementJacobian_u_u,
					 int eN)
    {
      for (int i=0;i<nDOF_test_element;i++)
	for (int j=0;j<nDOF_trial_element;j++)
	  {
	    elementJacobian_u_u[i*nDOF_trial_element+j]=0.0;
	  }
      double epsHeaviside,epsDirac,epsDiffusion;
      for  (int k=0;k<nQuadraturePoints_element;k++)
	{
	  int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
	    eN_k_nSpace = eN_k*nSpace;
	    //eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point
	  
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
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
	    }
	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      //int eN_k_i=eN_k*nDOF_test_element+i;
	      //int eN_k_i_nSpace=eN_k_i*nSpace;
	      int i_nSpace=i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  //int eN_k_j=eN_k*nDOF_trial_element+j;
		  //int eN_k_j_nSpace = eN_k_j*nSpace;
		  int j_nSpace = j*nSpace;
		  elementJacobian_u_u[i*nDOF_trial_element+j] += ck.ReactionJacobian_weak(1.0,
											  u_trial_ref[k*nDOF_trial_element+j],
											  u_test_dV[i]);
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
			   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			   double* globalJacobian)
    {
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double  elementJacobian_u_u[nDOF_test_element*nDOF_trial_element];
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      register int eN_j = eN*nDOF_trial_element+j;
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
				   elementJacobian_u_u,
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
    }//computeJacobian
  };//cppPres

  inline cppPres_base* newPres(int nSpaceIn,
                                 int nQuadraturePoints_elementIn,
                                 int nDOF_mesh_trial_elementIn,
                                 int nDOF_trial_elementIn,
                                 int nDOF_test_elementIn,
                                 int nQuadraturePoints_elementBoundaryIn,
                                 int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<cppPres_base,cppPres,CompKernel>(nSpaceIn,
                                                                                               nQuadraturePoints_elementIn,
                                                                                               nDOF_mesh_trial_elementIn,
                                                                                               nDOF_trial_elementIn,
                                                                                               nDOF_test_elementIn,
                                                                                               nQuadraturePoints_elementBoundaryIn,
                                                                                               CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<cppPres_base,cppPres,CompKernel>(nSpaceIn,
                                                                                             nQuadraturePoints_elementIn,
                                                                                             nDOF_mesh_trial_elementIn,
                                                                                             nDOF_trial_elementIn,
                                                                                             nDOF_test_elementIn,
                                                                                             nQuadraturePoints_elementBoundaryIn,
                                                                                             CompKernelFlag);
  }
}//proteus
#endif
