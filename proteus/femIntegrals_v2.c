#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include "femIntegrals_v2.h"
#include PROTEUS_LAPACK_H

/**
   \file femIntegrals.h
   \ingroup femIntegrals
   
   \brief A library of functions for computing the discrete finite element
   formulations.

@{
*/

void parametricFiniteElementSpace_getHessianValues(int nElements_global,
						   int nQuadraturePoints_element,
						   int nDOF_element,
						   int nSpace_global,
						   double* Hessian_psi,
						   double* inverseJacobianArray,
						   double* Hessian_vArray)
{
  int eN,k,j,I,J,K,L,nSpace_global2 = nSpace_global*nSpace_global;
  for(eN=0;eN<nElements_global;eN++)
    for(k=0;k<nQuadraturePoints_element;k++)
      for(j=0;j<nDOF_element;j++)
	for (I=0;I<nSpace_global;I++)
	  for (J=0;J<nSpace_global;J++)
	    for (K=0;K<nSpace_global;K++)
	      for (L=0;L<nSpace_global;L++)
		Hessian_vArray[eN*nQuadraturePoints_element*nDOF_element*nSpace_global2+
			       k*nDOF_element*nSpace_global2+
			       j*nSpace_global2+
			       I*nSpace_global+
			       J]
		  +=
		  Hessian_psi[k*nDOF_element*nSpace_global2+
			      j*nSpace_global2+
			      K*nSpace_global+
			      L]
		  *
		  inverseJacobianArray[eN*nQuadraturePoints_element*nSpace_global2+
				       k*nSpace_global2+
				       L*nSpace_global+
				       J]
		  *
		  inverseJacobianArray[eN*nQuadraturePoints_element*nSpace_global2+
				       k*nSpace_global2+
				       K*nSpace_global+
				       I];
}

void updateDiffusion2_strong(int nElements_global,
			     int nQuadraturePoints_element,
			     int nSpace,
			     double* a,
			     double* Hess_phi,
			     double* strong_residual)
{
  int eN,k,I,J,nSpace2=nSpace*nSpace;
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for(I=0;I<nSpace;I++)
	for (J=0;J<nSpace;J++)
	  {
            strong_residual[eN*nQuadraturePoints_element+
                            k] 
              -=
              a[eN*nQuadraturePoints_element*nSpace2 + 
                k*nSpace2 + 
                I*nSpace+
                J]
              *
              Hess_phi[eN*nQuadraturePoints_element*nSpace2 + 
                       k*nSpace2 +
                       J*nSpace+
                       I];
          }
}

void updateDiffusionJacobian2_strong(int nElements_global,
				     int nQuadraturePoints_element,
				     int nDOF_trial_element,
				     int nSpace,
				     int* l2g,
				     double* a,
				     double* da,
				     double* v,
				     double* Hess_phi,
				     double* dphi,
				     double* Hess_v,
				     double* dstrong_residual)
{
  int eN,k,j,I,J,nSpace2=nSpace*nSpace;
  /*double tmp;*/
  for(eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (j=0;j<nDOF_trial_element;j++)
	for(I=0;I<nSpace;I++)
	  for (J=0;J<nSpace;J++)
	    {
	      dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
			       k*nDOF_trial_element+
			       j]
		-=
		(da[eN*nQuadraturePoints_element*nSpace2 +
		    k*nSpace2 +
		    I*nSpace+
		    J]
		 *
		 v[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element+
		   j]
		 *
		 Hess_phi[eN*nQuadraturePoints_element*nSpace2 +
			  k*nSpace2 +
			  J*nSpace+
			  I]
		 +
		 a[eN*nQuadraturePoints_element*nSpace2 +
		   k*nSpace2 +
		   I*nSpace+
		   J]
		 *
		 dphi[l2g[eN*nDOF_trial_element +
			  j]]
		 *
		 Hess_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace2 +
			k*nDOF_trial_element*nSpace2 +
			j*nSpace2+
			J*nSpace+
			I]);
	    }
}

void updateDiffusion2_adjoint(int nElements_global,
			      int nQuadraturePoints_element,
			      int nDOF_test_element,
			      int nSpace,
			      double* a,
			      double* Hess_w_dV,
			      double* Lstar_w_dV)
{
  int eN,i,k,I,J,nSpace2=nSpace*nSpace;
  for(eN=0;eN<nElements_global;eN++)
    for (i=0;i<nDOF_test_element;i++)
      for (k=0;k<nQuadraturePoints_element;k++)
	for (I=0;I<nSpace;I++)
	  for(J=0;J<nSpace;J++)
            {
              Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
                         k*nDOF_test_element + 
                         i] 
                -=  
                a[eN*nQuadraturePoints_element*nSpace2 + 
                  k*nSpace2 + 
                  I*nSpace+
                  J]
                *
                Hess_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace2 + 
                          k*nDOF_test_element*nSpace2 +
                          i*nSpace2+
                          I*nSpace + 
                          J];
            }
}

void calculateWeightedShapeHessians(int nElements_global,
                                     int nQuadraturePoints_element,
                                     int nDOF_test_element,
                                     int nSpace,
                                     double* dVR,
                                     double* abs_det_jac,
                                     double* Hess_w,
                                     double* Hess_w_dV)
{
  int eN,i,k,I,J,nSpace2=nSpace*nSpace;
  for (eN=0;eN<nElements_global;eN++)
    for (k=0;k<nQuadraturePoints_element;k++)
      for (i=0;i<nDOF_test_element;i++)
        for (I=0;I<nSpace;I++)
	  for (J=0;J<nSpace;J++)
	    Hess_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace2+
		      k*nDOF_test_element*nSpace2+
		      i*nSpace2+
		      I*nSpace+
		      J] 
	      =
	      Hess_w[eN*nQuadraturePoints_element*nDOF_test_element*nSpace2+
		     k*nDOF_test_element*nSpace2+
		     i*nSpace2+
		     I*nSpace+
		     J]
	      *
	      dVR[k]
	      *
	      abs_det_jac[eN*nQuadraturePoints_element+
			  k];
}

void calculateFiniteElementFunctionHessianValues(int nElements_global,
						 int nQuadraturePoints_element,
						 int nDOF_trial_element,
						 int nComponents,
						 int nSpace,
						 int* l2g,
						 double* dof,
						 double* Hessian_v,
						 double* Hessian_u)
{
  int eN,k,j,t,I,J,nSpace2=nSpace*nSpace;
  memset(Hessian_u,0,sizeof(double)*nElements_global*nQuadraturePoints_element*nComponents*nSpace2);
  for (eN=0;eN<nElements_global;eN++)
    for  (k=0;k<nQuadraturePoints_element;k++)
      for (j=0;j<nDOF_trial_element;j++)
        for (t=0;t<nComponents;t++)
          for  (I=0;I<nSpace;I++)
	    for  (J=0;J<nSpace;J++)
	      Hessian_u[eN*nQuadraturePoints_element*nComponents*nSpace2+
			k*nComponents*nSpace2+
			t*nSpace2+
			I*nSpace+
			J]
		+=
		dof[l2g[eN*nDOF_trial_element+
			j]*nComponents+
		    t]
		*
		Hessian_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace2+
			  k*nDOF_trial_element*nSpace2+
			  j*nSpace2+
			  I*nSpace+
			  J];
}

/**
   \brief Update the global CSR Jacobian from the element boundary two-sided Hamiltonian flux Jacobians at interior boundaries
*/    
void updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_CSR(int nInteriorElementBoundaries_global,
									    int nElementBoundaries_element,
									    int nQuadraturePoints_elementBoundary,
									    int nDOF_test_element,
									    int nDOF_trial_element,
									    int* interiorElementBoundaries,
									    int* elementBoundaryElements,
									    int* elementBoundaryLocalElementBoundaries,
									    int* nFreeDOF_element_r,
									    int* freeLocal_r,
									    int* nFreeDOF_element_u,
									    int* freeLocal_u,
									    int* csrRowIndeces_ru,
									    int* csrColumnOffsets_eb_ru,
									    double* elementBoundaryFluxJacobian_2sided,
									    double* w_dS,
									    double* jac)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,ii,i,k,jj,j,nDOF_test_X_trial_element=nDOF_test_element*nDOF_trial_element,jacIndex;
  for (ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global   = elementBoundaryElements[ebN*2+0];
      right_eN_global  = elementBoundaryElements[ebN*2+1];
      left_ebN_element  = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      /*left flux*/
      for(ii=0;ii<nFreeDOF_element_r[left_eN_global];ii++)
        {
          i = freeLocal_r[left_eN_global*nDOF_test_element+ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element_u[left_eN_global];jj++)
                {
                  j = freeLocal_u[left_eN_global*nDOF_trial_element+
                                  jj];
                  jacIndex = csrRowIndeces_ru[left_eN_global*nDOF_test_element+
                                              ii]
                    +
                    csrColumnOffsets_eb_ru[ebN*4*nDOF_test_X_trial_element +
                                           0*2*nDOF_test_X_trial_element +
                                           0*nDOF_test_X_trial_element +
                                           ii*nDOF_trial_element +
                                           jj];
                  jac[jacIndex]
                    +=
                    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element +
						       0*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*left flux*/
						       0*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*left neig dep*/
						       k*nDOF_trial_element+
						       j]
                    *
                    w_dS[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element +
			 left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+
			 i];
                }
              for(jj=0;jj<nFreeDOF_element_u[right_eN_global];jj++)
                {
                  j = freeLocal_u[right_eN_global*nDOF_trial_element+
                                jj];
                  jacIndex = csrRowIndeces_ru[left_eN_global*nDOF_test_element+
                                           ii]
                    +
                    csrColumnOffsets_eb_ru[ebN*4*nDOF_test_X_trial_element+
                                        0*2*nDOF_test_X_trial_element+
                                        1*nDOF_test_X_trial_element+
                                        ii*nDOF_trial_element+
                                        jj];
                  jac[jacIndex]
                    +=
                    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element +
						       0*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*left flux*/
						       1*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*right neig dep*/
						       k*nDOF_trial_element+
						       j]
                    *
                    w_dS[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element +
			 left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+
			 i];
                }
            }
        }
      /*right flux*/
      for(ii=0;ii<nFreeDOF_element_r[right_eN_global];ii++)
        {
          i = freeLocal_r[right_eN_global*nDOF_test_element+
                        ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element_u[left_eN_global];jj++)
                {
                  j = freeLocal_u[left_eN_global*nDOF_trial_element+
                                jj];
                  jacIndex = csrRowIndeces_ru[right_eN_global*nDOF_test_element+
                                           ii]
                    +
                    csrColumnOffsets_eb_ru[ebN*4*nDOF_test_X_trial_element+
                                        1*2*nDOF_test_X_trial_element+
                                        0*nDOF_test_X_trial_element+
                                        ii*nDOF_trial_element+
                                        jj];
                  jac[jacIndex]
                    +=
                    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element+
						       1*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*right flux*/
						       0*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*left neig dep*/
						       k*nDOF_trial_element+
						       j]
                    *
                    w_dS[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+
			 i];
                }
              for(jj=0;jj<nFreeDOF_element_u[right_eN_global];jj++)
                {
                  j = freeLocal_u[right_eN_global*nDOF_trial_element+
                                jj];
                  jacIndex = csrRowIndeces_ru[right_eN_global*nDOF_test_element+
                                           ii]
                    +
                    csrColumnOffsets_eb_ru[ebN*4*nDOF_test_X_trial_element+
                                        1*2*nDOF_test_X_trial_element+
                                        1*nDOF_test_X_trial_element+
                                        ii*nDOF_trial_element+
                                        jj];
                  jac[jacIndex]
                    +=
                    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element+
						       1*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*right flux*/
						       1*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*right neig dep*/
						       k*nDOF_trial_element+
						       j]*
                    w_dS[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+
			 i];
                }
            }
        }
    }
}

/**
   \brief Update the global dense Jacobian from the element boundary two-sided Hamiltonflux Jacobians at interior boundaries
*/    
void updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_dense(int nInteriorElementBoundaries_global,
									  int nElementBoundaries_element,
									  int nQuadraturePoints_elementBoundary,
									  int nDOF_test_element,
									  int nDOF_trial_element,
									  int offset_r,
									  int stride_r,
									  int offset_u,
									  int stride_u,
									  int nFreeVDOF_global,
									  int* interiorElementBoundaries,
									  int* elementBoundaryElements,
									  int* elementBoundaryLocalElementBoundaries,
									  int* nFreeDOF_element_r,
									  int* nFreeDOF_element_u,
									  int* freeLocal_r,
									  int* freeGlobal_r,
									  int* freeLocal_u,
									  int* freeGlobal_u,
									  double* elementBoundaryFluxJacobian_2sided,
									  double* w_dS,
									  double* jac)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,ii,i,k,jj,j,jacIndex,I,J;
  for (ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global   = elementBoundaryElements[ebN*2+0];
      right_eN_global  = elementBoundaryElements[ebN*2+1];
      left_ebN_element  = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      /*left  flux*/
      for(ii=0;ii<nFreeDOF_element_r[left_eN_global];ii++)
	{
	  i = freeLocal_r[left_eN_global*nDOF_test_element+
			  ii];
	  I = offset_r + stride_r*freeGlobal_r[left_eN_global*nDOF_test_element+
					       ii];
	  for(k=0;k<nQuadraturePoints_elementBoundary;k++)
	    {
	      for(jj=0;jj<nFreeDOF_element_u[left_eN_global];jj++)
		{
		  j = freeLocal_u[left_eN_global*nDOF_trial_element+
				  jj];
		  J = offset_u + stride_u*freeGlobal_u[left_eN_global*nDOF_trial_element+
						       jj];
		  jacIndex = I+J*nFreeVDOF_global;
		  jac[jacIndex] 
		    += 
		    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
						       0*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*left flux*/
						       0*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*left  neig. dep*/
						       k*nDOF_trial_element+
						       j]
		    *
		    w_dS[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element + 
			 left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+
			 i];
		}
	      for(jj=0;jj<nFreeDOF_element_u[right_eN_global];jj++)
		{
		  j = freeLocal_u[right_eN_global*nDOF_trial_element+
				  jj];
		  J = offset_u + stride_u*freeGlobal_u[right_eN_global*nDOF_trial_element+
						       jj];
		  jacIndex = I+J*nFreeVDOF_global;
		  jac[jacIndex] 
		    += 
		    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + 
						       0*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*left flux*/
						       1*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*right  neig. dep*/
						       k*nDOF_trial_element+
						       j]*
		    w_dS[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element + 
			 left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+i];
		}
	    }/*k*/
	}/*ii*/
      /*right flux*/
      for(ii=0;ii<nFreeDOF_element_r[right_eN_global];ii++)
        {
          i = freeLocal_r[right_eN_global*nDOF_test_element+
                        ii];
	  I = offset_r + stride_r*freeGlobal_r[right_eN_global*nDOF_test_element+
					     ii];
          for(k=0;k<nQuadraturePoints_elementBoundary;k++)
            {
              for(jj=0;jj<nFreeDOF_element_u[left_eN_global];jj++)
                {
                  j = freeLocal_u[left_eN_global*nDOF_trial_element+
                                jj];
		  J = offset_u + stride_u*freeGlobal_u[left_eN_global*nDOF_trial_element+
						     jj];
                  jacIndex = I+J*nFreeVDOF_global;
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element+
						       1*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*right flux*/
						       0*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*left neighbor*/
						       k*nDOF_trial_element+
						       j]
                    *
                    w_dS[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+
			 i];
                }
              for(jj=0;jj<nFreeDOF_element_u[right_eN_global];jj++)
                {
                  j = freeLocal_u[right_eN_global*nDOF_trial_element+
                                jj];
		  J = offset_u + stride_u*freeGlobal_u[right_eN_global*nDOF_trial_element+
						     jj];
                  jacIndex = I+J*nFreeVDOF_global;
                  jac[jacIndex] 
                    += 
                    elementBoundaryFluxJacobian_2sided[ebN*2*2*nQuadraturePoints_elementBoundary*nDOF_trial_element+
						       1*2*nQuadraturePoints_elementBoundary*nDOF_trial_element + /*right flux*/
						       1*nQuadraturePoints_elementBoundary*nDOF_trial_element+ /*right neig dep*/
						       k*nDOF_trial_element+
						       j]*
                    w_dS[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
			 k*nDOF_test_element+
			 i];
                }
            }
        }
    }
}

/**
   \brief Update a two-sided (say nonconservative HJ flux) element boundary flux on interior element boundaries
*/
void updateInteriorTwoSidedElementBoundaryFlux(int nInteriorElementBoundaries_global,
					       int nElementBoundaries_element,
					       int nQuadraturePoints_elementBoundary,
					       int nDOF_test_element,
					       int* interiorElementBoundaries,
					       int* elementBoundaryElements,
					       int* elementBoundaryLocalElementBoundaries,
					       double* flux,
					       double* w_dS,
					       double* residual)
{
  int ebNI,ebN,left_eN_global,right_eN_global,left_ebN_element,right_ebN_element,i,k;
  for(ebNI=0;ebNI<nInteriorElementBoundaries_global;ebNI++)
    {
      ebN = interiorElementBoundaries[ebNI];
      left_eN_global = elementBoundaryElements[ebN*2+0];
      right_eN_global = elementBoundaryElements[ebN*2+1];
      left_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      right_ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+1];
      for(i=0;i<nDOF_test_element;i++)
        for(k=0;k<nQuadraturePoints_elementBoundary;k++)
          {
            residual[left_eN_global*nDOF_test_element+
		     i]
              +=
              flux[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
		   left_ebN_element*nQuadraturePoints_elementBoundary+ /*left flux*/
                   k]*
              w_dS[left_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
		   left_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
		   k*nDOF_test_element+
		   i];
            residual[right_eN_global*nDOF_test_element+
		     i]
              +=
              flux[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
		   right_ebN_element*nQuadraturePoints_elementBoundary+ /*right flux*/
                   k]*
              w_dS[right_eN_global*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
		   right_ebN_element*nQuadraturePoints_elementBoundary*nDOF_test_element+
		   k*nDOF_test_element+
		   i];
          }
    }
}

/*not really likely but in case have two separate characteristic speeds to add?*/
void calculateCFLADR2speeds(int nElements_global,
			    int nQuadraturePoints_element,
			    int nSpace,
			    double* elementDiameter,
			    double* dm,
			    double* df1,
			    double* df2,
			    double* cfl)
{
  int eN,k,I;
  double h,Vlin,dmdu,cfl1;
  for(eN=0;eN<nElements_global;eN++)
    {
      h = elementDiameter[eN];
      for (k=0;k<nQuadraturePoints_element;k++)
        {
	  dmdu = dm[eN*nQuadraturePoints_element + k];
          Vlin = 0.0;
          for(I=0;I<nSpace;I++)
            Vlin 
              += 
              df1[eN*nQuadraturePoints_element*nSpace + 
                 k*nSpace + 
                 I]
              *
              df1[eN*nQuadraturePoints_element*nSpace + 
		  k*nSpace + 
		  I]
	      + 
              df2[eN*nQuadraturePoints_element*nSpace + 
		  k*nSpace + 
		  I]
              *
              df2[eN*nQuadraturePoints_element*nSpace + 
		  k*nSpace + 
		  I];

          Vlin = sqrt(Vlin)/(dmdu+1.0e-10);
          cfl1 = Vlin/h;
	  cfl[eN*nQuadraturePoints_element + 
	      k] = cfl1;

        }
    }
}

/**
  for debugging
 */
int checkElementBoundaryAndExteriorElementBoundaryArraysSame(int nElementBoundaries_element,
							     int nExteriorElementBoundaries_global,
							     int nQuadraturePoints_elementBoundary,
							     int nValuesPerQuadraturePoint,
							     double tolerance,
							     const int * exteriorElementBoundariesArray,
							     const int * elementBoundaryElementsArray,
							     const int * elementBoundaryLocalElementBoundariesArray,
							     const double * ebq_val,
							     const double * ebqe_val,
							     int* firstBadIndex)
{
  int eN,ebN,ebN_local,ebNE,k,i;
  double val1,val2;
  int failed = 0;
  *firstBadIndex = -1;

  for (ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      ebN = exteriorElementBoundariesArray[ebNE];
      eN  = elementBoundaryElementsArray[ebN*2+0];
      ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0];

      for (i=0; i < nValuesPerQuadraturePoint; i++)
	{
	  val1 = ebq_val[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nValuesPerQuadraturePoint +
			 ebN_local*nQuadraturePoints_elementBoundary*nValuesPerQuadraturePoint +
			 k*nValuesPerQuadraturePoint + 
			 i];
	  val2 = ebqe_val[ebNE*nQuadraturePoints_elementBoundary*nValuesPerQuadraturePoint +
			  k*nValuesPerQuadraturePoint + 
			  i];
	  if (fabs(val1-val2) > tolerance)
	    {
	      failed = 1;
	      *firstBadIndex = ebNE;
	      return failed;
	    }
	}
    }
  return failed;
}

/**
  for debugging
 */
int checkGlobalElementBoundaryAndExteriorElementBoundaryArraysSame(int nExteriorElementBoundaries_global,
								   int nQuadraturePoints_elementBoundary,
								   int nValuesPerQuadraturePoint,
								   double tolerance,
								   const int * exteriorElementBoundariesArray,
								   const int * elementBoundaryElementsArray,
								   const int * elementBoundaryLocalElementBoundariesArray,
								   const double * ebq_global_val,
								   const double * ebqe_val,
								   int* firstBadIndex)
{
  int eN,ebN,ebNE,k,i;
  double val1,val2;
  int failed = 0;
  *firstBadIndex = -1;

  for (ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      ebN = exteriorElementBoundariesArray[ebNE];
      eN  = elementBoundaryElementsArray[ebN*2+0];

      for (i=0; i < nValuesPerQuadraturePoint; i++)
	{
	  val1 = ebq_global_val[ebN*nQuadraturePoints_elementBoundary*nValuesPerQuadraturePoint +
				k*nValuesPerQuadraturePoint + 
				i];
	  val2 = ebqe_val[ebNE*nQuadraturePoints_elementBoundary*nValuesPerQuadraturePoint +
			  k*nValuesPerQuadraturePoint + 
			  i];
	  if (fabs(val1-val2) > tolerance)
	    {
	      failed = 1;
	      *firstBadIndex = ebNE;
	      return failed;
	    }
	}
    }
  return failed;
}

void calculateExteriorElementBoundaryStress3D(int nExteriorElementBoundaries_global,
                                              int nQuadraturePoints_elementBoundary,
                                              int* elementBoundaryMaterialTypes,
                                              int* exteriorElementBoundaries,
                                              int* elementBoundaryElements,
                                              int* elementBoundaryLocalElementBoundaries,
                                              double* p,
                                              double* mom_flux_vec_u,
                                              double* mom_flux_vec_v,
                                              double* mom_flux_vec_w,
                                              double* dS,
                                              double* n,
                                              double* F)
{
  int ebNE,ebN,eN_global,ebN_element,k,I;
  double stress_x,stress_y,stress_z;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          stress_x = 0.0;
          stress_y = 0.0;
          stress_z = 0.0;
          stress_x += p[ebNE*nQuadraturePoints_elementBoundary+
                        k]
            *
            n[ebNE*nQuadraturePoints_elementBoundary*3+
              k*3+
              0];
          stress_y += p[ebNE*nQuadraturePoints_elementBoundary+
                        k]
            *
            n[ebNE*nQuadraturePoints_elementBoundary*3+
              k*3+
              1];
          stress_z += p[ebNE*nQuadraturePoints_elementBoundary+
                        k]
            *
            n[ebNE*nQuadraturePoints_elementBoundary*3+
              k*3+
              2];
          for(I=0;I<3;I++)
            {
              stress_x += mom_flux_vec_u[ebNE*nQuadraturePoints_elementBoundary*3+
                                         k*3+
                                         I]
                *
                n[ebNE*nQuadraturePoints_elementBoundary*3+
                  k*3+
                  I];
              stress_y += mom_flux_vec_v[ebNE*nQuadraturePoints_elementBoundary*3+
                                         k*3+
                                         I]
                *
                n[ebNE*nQuadraturePoints_elementBoundary*3+
                  k*3+
                  I];
              stress_z += mom_flux_vec_w[ebNE*nQuadraturePoints_elementBoundary*3+
                                         k*3+
                                         I]
                *
                n[ebNE*nQuadraturePoints_elementBoundary*3+
                  k*3+
                  I];
            }
          F[elementBoundaryMaterialTypes[ebN]*3 + 0] += stress_x*
            dS[ebNE*nQuadraturePoints_elementBoundary+
               k];
          F[elementBoundaryMaterialTypes[ebN]*3 + 1] += stress_y*
            dS[ebNE*nQuadraturePoints_elementBoundary+
               k];
          F[elementBoundaryMaterialTypes[ebN]*3 + 2] += stress_z*
            dS[ebNE*nQuadraturePoints_elementBoundary+
               k];
        }
    }
}

void calculateExteriorElementBoundaryStress2D(int nExteriorElementBoundaries_global,
                                              int nQuadraturePoints_elementBoundary,
                                              int* elementBoundaryMaterialTypes,
                                              int* exteriorElementBoundaries,
                                              int* elementBoundaryElements,
                                              int* elementBoundaryLocalElementBoundaries,
                                              double* p,
                                              double* mom_flux_vec_u,
                                              double* mom_flux_vec_v,
                                              double* dS,
                                              double* n,
                                              double* F)
{
  int ebNE,ebN,eN_global,ebN_element,k,I;
  double stress_x,stress_y;
  for(ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      ebN = exteriorElementBoundaries[ebNE];
      eN_global = elementBoundaryElements[ebN*2+0];
      ebN_element = elementBoundaryLocalElementBoundaries[ebN*2+0];
      for(k=0;k<nQuadraturePoints_elementBoundary;k++)
        {
          stress_x = 0.0;
          stress_y = 0.0;
          stress_x += p[ebNE*nQuadraturePoints_elementBoundary+
                        k]
            *
            n[ebNE*nQuadraturePoints_elementBoundary*2+
              k*2+
              0];
          stress_y += p[ebNE*nQuadraturePoints_elementBoundary+
                        k]
            *
            n[ebNE*nQuadraturePoints_elementBoundary*2+
              k*2+
              1];
          /* for(I=0;I<2;I++) */
          /*   { */
          /*     stress_x += mom_flux_vec_u[ebNE*nQuadraturePoints_elementBoundary*2+ */
          /*                                k*2+ */
          /*                                I] */
          /*       * */
          /*       n[ebNE*nQuadraturePoints_elementBoundary*2+ */
          /*         k*2+ */
          /*         I]; */
          /*     stress_y += mom_flux_vec_v[ebNE*nQuadraturePoints_elementBoundary*2+ */
          /*                                k*2+ */
          /*                                I] */
          /*       * */
          /*       n[ebNE*nQuadraturePoints_elementBoundary*2+ */
          /*         k*2+ */
          /*         I]; */
          /*   } */
          F[elementBoundaryMaterialTypes[ebN]*2 + 0] += stress_x*
            dS[ebNE*nQuadraturePoints_elementBoundary+
               k];
          F[elementBoundaryMaterialTypes[ebN]*2 + 1] += stress_y*
            dS[ebNE*nQuadraturePoints_elementBoundary+
               k];
        }
    }
}
