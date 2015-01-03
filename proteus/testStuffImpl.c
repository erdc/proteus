#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <assert.h>
#include "testStuffImpl.h"
#include PROTEUS_LAPACK_H

/** \file testStuffImpl.c
    \ingroup testStuffImpl
    @{
*/

/***********************************************************************
  Apply Algorithm 5.2 in Barth and Sethian 98 for level set equation

  nElements_global : IN, number of elements (Ne)
  nDOF_element     : IN, number of dofs per element := nSpace+1
  nSpace           : IN, space dimension (Nsd)
  nQuadraturePoints_element: IN, number of quadrature points per element (nq)
  nQuadraturePoints_elementBoundary : IN, number of quadrature points per
                    element boundary (nqb)

  elementDiameters : IN, Ne x 1. 
                     elementDiameters[I]=h_{E_I} for global element I
  elementUnitOuterNormals : IN, Ne x Nsd+1 x nqb x Nsd. unit outer normals for element boundaries
                     elementUnitOuterNormals[I,j,k,:] = \vec n^j_I(x^j_k)
  detJ             : IN, Ne x nq.  
                     detJ[I,k] = det J_{I,k} for global element I and 
                     quadrature point k
                     = |E_I|/[1,2,6] for 1d,2d,3d
  sqrtDetG         : IN, Ne x Nsd+1 x nqb. sqrt of metric tensor deterimant
                     for local boundary j of global element I eval'd at 
                     quadrature point k
                     sqrtDetG[I,j,k] = sqrt(det(g_I^j)(x_k))  
                     = |e_I^j|/[1,1,2]  for 1d,2d,3d
  elementQuadratureWeights : IN, nq x 1. 
                     element quadrature points on reference element

  l2g              : IN, Ne x ndof. local to global degree of freedom map
                     l2g[I,j] = global dof for element I's j'th degree of freedom
  phiIn            : IN, Ndof x 1. input degrees of freeom.
 
  H                : IN, Ne x nq. Hamiltonian F|\grad phi| eval'd at element
                     quad points
                     H[I,k] = H(\vec x^k_I)

  dH               : IN, Ne x nq x Nsd. derivative of H w.r.t \grad phi.
                     dH[I,k,:]=\pd{H}{\grad \phi}(\vec x^k_{I})

  r                : IN, Ne x nq. minus source term. 
                     r[I,k]  = -f(x^k_{I})
  dt               : IN, \Delta t
  phiOut           : IN, Ndof x 1. output degrees of freeom.
  wOut             : IN, Ndof x 1. output test function weights
  
 ***********************************************************************/


void advanceStageP1_C0_GLS_lump(int nElements_global,
				int nElementBoundaries_element,
				int nDOF_element,
				int nSpace,
				int nQuadraturePoints_element,
				int nQuadraturePoints_elementBoundary,
				double* elementDiameters,
				double* elementUnitOuterNormals,
				double* detJ,
				double* sqrtDetG,
				double* elementQuadratureWeights,
				int* l2g,
				double * phiIn,
				double * H,
				double *dH,
				double *r,
				double dt,
				double *phiOut,
				double *wOut)
{
  int eN,i,k,id;
  double vol_eN,area_i,h_eN,d,volFact,areaFact,ntmp;
  double Hbar,fbar,gphiNorm,gradHbar,tau,deltaPhi,ai,nid;
  double gradPhi[3] = {0.,0.,0.};
  double K[4] = {0.,0.,0.,0.};
  assert(nElementBoundaries_element == nSpace+1);
  assert(nDOF_element == nSpace+1);
  d = nSpace;
  /*mwf debug 
  printf("Entering advance LS eN=%d neb=%d nde=%d Nsd=%d nq=%d nqb=%d\n",
	 nElements_global,nElementBoundaries_element,nDOF_element,
	 nSpace,nQuadraturePoints_element,nQuadraturePoints_elementBoundary);
  */
  /*convert metric tensors to volumes and areas for formulas*/
  volFact = 1.0; areaFact = 1.0;
  if (nSpace > 1)
    volFact = 0.5;
  if (nSpace > 2)
    { volFact = 1.0/6.0; areaFact = 0.5;}

  /*zero phiOut first*/
  for (eN = 0; eN < nElements_global; eN++)
    for (i=0; i < nDOF_element; i++)
      {
	phiOut[l2g[eN*nDOF_element+i]] = 0.0;
	wOut[l2g[eN*nDOF_element+i]] = 0.0;
      }
  for (eN = 0; eN < nElements_global; eN++)
    {
      h_eN = elementDiameters[eN];
      vol_eN = fabs(detJ[eN*nQuadraturePoints_element + 0])*volFact;/*assumed affine*/
      
      /*element averages using quantities already computed at quad points*/
      Hbar = 0.0; fbar = 0.0; gradHbar = 0.0;
      for (k = 0; k < nQuadraturePoints_element; k++)
	{
	  Hbar += H[eN*nQuadraturePoints_element+k]*elementQuadratureWeights[k]*
	    vol_eN/volFact;
	  fbar += -r[eN*nQuadraturePoints_element+k]*elementQuadratureWeights[k]*
	    vol_eN/volFact;
	  ntmp = 0.0;
	  for (id = 0; id < nSpace; id++)
	    {
	      ntmp  += 
		dH[eN*nQuadraturePoints_element*nSpace + 
		   k*nSpace + id]
		 *
		 dH[eN*nQuadraturePoints_element*nSpace + 
		    k*nSpace + id];
	      /*mwf debug
	      printf("dH[%d,%d,%d]=%g \n",eN,k,id,
		     dH[eN*nQuadraturePoints_element*nSpace +k*nSpace + id]);
	      */
	    }
	  gradHbar += sqrt(ntmp)*elementQuadratureWeights[k]*vol_eN/volFact;
	}/*k*/
      Hbar = Hbar/vol_eN; fbar = fbar/vol_eN; gradHbar = gradHbar/vol_eN;
      /*mwf debug
      printf("eN=%d vol=%g Hbar=%g fbar=%g gradHbar=%g\n",eN,vol_eN,
	     Hbar,fbar,gradHbar);
      */
      /*compute gradient of phiIn on element (constant)*/
      gradPhi[0] = 0.; gradPhi[1] = 0.; gradPhi[2] = 0.;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      /*nid is inner normal with area*/
	      nid = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradPhi[id] += phiIn[l2g[eN*nDOF_element+i]]*nid/(d*vol_eN);
	      /*mwf debug
	      printf("eN=%d ei=%i area=%g nid=%g gradPhi[%d]=%g\n",eN,i,area_i,nid,id,gradPhi[id]);
	      */
	    }/*end id*/
	  
	}/*i*/
      gphiNorm = 0.0;
      for (id = 0; id < nSpace; id++)
	gphiNorm += gradPhi[id]*gradPhi[id];
      gphiNorm = sqrt(gphiNorm);
      /*mwf debug
      printf("eN=%d vol=%g Hbar=%g gradHbar=%g gradPhi[0]=%g  Hbar/|gradPhi|= %g \n",eN,vol_eN,
	     Hbar,gradHbar,gradPhi[0],Hbar/(gphiNorm+1.0e-12));
      */
      
      /*tau formula from eqn (58) */
      tau = sqrt(4.0/(dt*dt) + 4.0*gradHbar*gradHbar/(h_eN*h_eN) + 1.0e-12);
      /*now try with 1/dt instead of 2/dt*/
      /*tau = sqrt(1.0/(dt*dt) + 4.0*gradHbar*gradHbar/(h_eN*h_eN) + 1.0e-12);*/
      tau = 1.0/tau;

      deltaPhi = 0.0;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  K[i] = 0.0; /*save for ai calc*/
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      nid = -area_i /*inner normal with area*/
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /*for LS Fbar = Hbar/|\grad phi| */
	      K[i] += Hbar*gradPhi[id]*nid/(d*gphiNorm*gphiNorm+1.0e-12);
	      /*mwf debug
	      printf("eN=%d i=%i Hbar=%g gradPhi[0]=%g K[%d]=%g \n",eN,i,Hbar,gradPhi[0],i,K[i]);
	      */
	    } /*id*/
	  deltaPhi += phiIn[l2g[eN*nDOF_element+i]]*K[i];
	  
	}/*i*/
      /*final update*/
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  ai = 1./d + tau*K[i];
	  phiOut[l2g[eN*nDOF_element+i]] += ai*(deltaPhi - fbar*vol_eN);
	  wOut[l2g[eN*nDOF_element+i]]   += ai*vol_eN;
	}
    }/*eN*/ 
  
}

void advanceStageP1_C0_GLS_lump_noSource(int nElements_global,
					 int nElementBoundaries_element,
					 int nDOF_element,
					 int nSpace,
					 int nQuadraturePoints_element,
					 int nQuadraturePoints_elementBoundary,
					 double* elementDiameters,
					 double* elementUnitOuterNormals,
					 double* detJ,
					 double* sqrtDetG,
					 double* elementQuadratureWeights,
					 int* l2g,
					 double * phiIn,
					 double * H,
					 double *dH,
					 double dt,
					 double *phiOut,
					 double *wOut)
{
  int eN,i,k,id;
  double vol_eN,area_i,h_eN,d,volFact,areaFact,ntmp;
  double Hbar,gphiNorm,gradHbar,tau,deltaPhi,ai,nid;
  double gradPhi[3] = {0.,0.,0.};
  double K[4] = {0.,0.,0.,0.};
  assert(nElementBoundaries_element == nSpace+1);
  assert(nDOF_element == nSpace+1);
  d = nSpace;
  /*mwf debug 
  printf("Entering advance LS eN=%d neb=%d nde=%d Nsd=%d nq=%d nqb=%d\n",
	 nElements_global,nElementBoundaries_element,nDOF_element,
	 nSpace,nQuadraturePoints_element,nQuadraturePoints_elementBoundary);
  */
  /*convert metric tensors to volumes and areas for formulas*/
  volFact = 1.0; areaFact = 1.0;
  if (nSpace > 1)
    volFact = 0.5;
  if (nSpace > 2)
    { volFact = 1.0/6.0; areaFact = 0.5;}

  /*zero phiOut first*/
  for (eN = 0; eN < nElements_global; eN++)
    for (i=0; i < nDOF_element; i++)
      {
	phiOut[l2g[eN*nDOF_element+i]] = 0.0;
	wOut[l2g[eN*nDOF_element+i]] = 0.0;
      }
  for (eN = 0; eN < nElements_global; eN++)
    {
      h_eN = elementDiameters[eN];
      vol_eN = fabs(detJ[eN*nQuadraturePoints_element + 0])*volFact;/*assumed affine*/
      
      /*element averages using quantities already computed at quad points*/
      Hbar = 0.0; gradHbar = 0.0;
      for (k = 0; k < nQuadraturePoints_element; k++)
	{
	  Hbar += H[eN*nQuadraturePoints_element+k]*elementQuadratureWeights[k]*
	    vol_eN/volFact;
	  ntmp = 0.0;
	  for (id = 0; id < nSpace; id++)
	    {
	      ntmp  += 
		dH[eN*nQuadraturePoints_element*nSpace + 
		   k*nSpace + id]
		 *
		 dH[eN*nQuadraturePoints_element*nSpace + 
		    k*nSpace + id];
	      /*mwf debug
	      printf("dH[%d,%d,%d]=%g \n",eN,k,id,
		     dH[eN*nQuadraturePoints_element*nSpace +k*nSpace + id]);
	      */
	    }
	  gradHbar += sqrt(ntmp)*elementQuadratureWeights[k]*vol_eN/volFact;
	}/*k*/
      Hbar = Hbar/vol_eN; gradHbar = gradHbar/vol_eN;
      /*mwf debug
      printf("eN=%d vol=%g Hbar=%g gradHbar=%g\n",eN,vol_eN,
	     Hbar,gradHbar);
      */
      /*compute gradient of phiIn on element (constant)*/
      gradPhi[0] = 0.; gradPhi[1] = 0.; gradPhi[2] = 0.;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      /*nid is inner normal with area*/
	      nid = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradPhi[id] += phiIn[l2g[eN*nDOF_element+i]]*nid/(d*vol_eN);
	      /*mwf debug
	      printf("eN=%d ei=%i area=%g nid=%g gradPhi[%d]=%g\n",eN,i,area_i,nid,id,gradPhi[id]);
	      */
	    }/*end id*/
	  
	}/*i*/
      gphiNorm = 0.0;
      for (id = 0; id < nSpace; id++)
	gphiNorm += gradPhi[id]*gradPhi[id];
      gphiNorm = sqrt(gphiNorm);
      /*mwf debug
      printf("eN=%d vol=%g Hbar=%g gradHbar=%g gradPhi[0]=%g  Hbar/|gradPhi|= %g \n",eN,vol_eN,
	     Hbar,gradHbar,gradPhi[0],Hbar/(gphiNorm+1.0e-12));
      */
      
      /*tau formula from eqn (58) */
      tau = sqrt(4.0/(dt*dt) + 4.0*gradHbar*gradHbar/(h_eN*h_eN) + 1.0e-12);
      /*now try with 1/dt instead of 2/dt*/
      /*tau = sqrt(1.0/(dt*dt) + 4.0*gradHbar*gradHbar/(h_eN*h_eN) + 1.0e-12)*/;
      /*now try with 1 norm*/
      /*tau = 1.0/dt + fabs(2.0*gradHbar/h_eN) + 1.0e-12;*/
      tau = 1.0/tau;
      /*mwf debug set tau to zero? didn't change much*/
      deltaPhi = 0.0;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  K[i] = 0.0; /*save for ai calc*/
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      nid = -area_i /*inner normal with area*/
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /*for LS Fbar = Hbar/|\grad phi| */
	      K[i] += Hbar*gradPhi[id]*nid/(d*gphiNorm*gphiNorm+1.0e-12);
	      /*mwf debug
	      printf("eN=%d i=%i Hbar=%g gradPhi[0]=%g K[%d]=%g \n",eN,i,Hbar,gradPhi[0],i,K[i]);
	      */
	    } /*id*/
	  deltaPhi += phiIn[l2g[eN*nDOF_element+i]]*K[i];
	  
	}/*i*/
      /*final update*/
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  /*mwf why not 1/(d+1)?*/
	  ai = 1./d + tau*K[i];
	  phiOut[l2g[eN*nDOF_element+i]] += ai*deltaPhi;
	  wOut[l2g[eN*nDOF_element+i]]   += ai*vol_eN;
	}
    }/*eN*/ 
  
}


void advanceStageP1_C0_SGS_lump_noSource(int nElements_global,
					 int nElementBoundaries_element,
					 int nDOF_element,
					 int nSpace,
					 int nQuadraturePoints_element,
					 int nQuadraturePoints_elementBoundary,
					 double* elementDiameters,
					 double* elementUnitOuterNormals,
					 double* detJ,
					 double* sqrtDetG,
					 double* elementQuadratureWeights,
					 int* l2g,
					 double * phiIn,
					 double * H,
					 double *dH,
					 double dt,
					 double *phiOut,
					 double *wOut)
{
  int eN,i,id;
  double vol_eN,area_i,h_eN,d,volFact,areaFact;
  double nid,gphiNorm,gradHnorm,gradHdotGradNi,gradHdotGradPhi,tau,
    dhid,ai,gradNid;
  double gradPhi[3] = {0.,0.,0.};
  double K[4] = {0.,0.,0.,0.};

  assert(nElementBoundaries_element == nSpace+1);
  assert(nDOF_element == nSpace+1);
  d = nSpace;
  /*mwf debug 
  printf("Entering advance LS SGS eN=%d neb=%d nde=%d Nsd=%d nq=%d nqb=%d\n",
	 nElements_global,nElementBoundaries_element,nDOF_element,
	 nSpace,nQuadraturePoints_element,nQuadraturePoints_elementBoundary);
  */
  /*convert metric tensors to volumes and areas for formulas*/
  volFact = 1.0; areaFact = 1.0;
  if (nSpace > 1)
    volFact = 0.5;
  if (nSpace > 2)
    { volFact = 1.0/6.0; areaFact = 0.5;}

  /*zero phiOut first*/
  for (eN = 0; eN < nElements_global; eN++)
    for (i=0; i < nDOF_element; i++)
      {
	phiOut[l2g[eN*nDOF_element+i]] = 0.0;
	wOut[l2g[eN*nDOF_element+i]] = 0.0;
      }
  for (eN = 0; eN < nElements_global; eN++)
    {
      h_eN = elementDiameters[eN];
      vol_eN = fabs(detJ[eN*nQuadraturePoints_element + 0])*volFact;/*assumed affine*/
      
      /*compute gradient of phiIn on element (constant)*/
      gradPhi[0] = 0.; gradPhi[1] = 0.; gradPhi[2] = 0.;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      /*nid is inner normal with area*/
	      nid = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradPhi[id] += phiIn[l2g[eN*nDOF_element+i]]*nid/(d*vol_eN);
	      /*mwf debug
	      printf("eN=%d ei=%i area=%g nid=%g gradPhi[%d]=%g\n",eN,i,area_i,nid,id,gradPhi[id]);
	      */
	    }/*end id*/
	  
	}/*i*/
      gphiNorm = 0.0;
      for (id = 0; id < nSpace; id++)
	gphiNorm += gradPhi[id]*gradPhi[id];
      gphiNorm = sqrt(gphiNorm);

      /*
	compute K[i] = \tau_i (\grad H_i \dot \grad N_i)(\grad H_i \dot \grad phi)
	relies on vertex quadrature used to eval dH and correspondence with loc dofs
      */
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/

	  gradHnorm = 0.0; gradHdotGradNi = 0.0; gradHdotGradPhi=0.0;
	  for (id = 0; id < nSpace; id++)
	    {
	      dhid = dH[eN*nQuadraturePoints_element*nSpace + i*nSpace + id];
	      nid  = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradNid = nid/(d*vol_eN);
	      
	      gradHnorm      += dhid*dhid;
	      gradHdotGradNi += dhid*gradNid;
	      gradHdotGradPhi+= dhid*gradPhi[id];

	    }
	  gradHnorm = sqrt(gradHnorm);
	  
	  /*tau formula from eqn (58) */
	  tau = sqrt(4.0/(dt*dt) + 4.0*gradHnorm*gradHnorm/(h_eN*h_eN) + 1.0e-12);
	  /*try 1 norm */
	  /*tau = 1./dt + 2.*gradHnorm/h_eN + 1.0e-12;*/
	  tau = 1.0/tau;

	  /*mwf debug
	  if (tau < 1.0e-4 || fabs(gradHdotGradNi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradNi=%g \n",eN,i,gradHdotGradNi);
	  if (fabs(gradHdotGradPhi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradPhi=%g \n",eN,i,gradHdotGradPhi);
	  */
	  /*mwf debug sharpens result for L_S_US_1d_circle 
	    tau = tau*0.01;
	  */
	  
	  K[i] = (1.0 + tau*gradHdotGradNi)*gradHdotGradPhi;
	}/* i */
      /*final update*/
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  ai = vol_eN/(d+1.0);
	  phiOut[l2g[eN*nDOF_element+i]] += ai*K[i]; 
	  wOut[l2g[eN*nDOF_element+i]]   += ai;
	}
    }/*eN*/ 
  
}

void advanceStageP1_C0_SUPG_lump(int nElements_global,
				 int nElementBoundaries_element,
				 int nDOF_element,
				 int nSpace,
				 int nQuadraturePoints_element,
				 int nQuadraturePoints_elementBoundary,
				 double* elementDiameters,
				 double* elementUnitOuterNormals,
				 double* detJ,
				 double* sqrtDetG,
				 double* elementQuadratureWeights,
				 int* l2g,
				 double * phiIn,
				 double * H,
				 double *dH,
				 double *r,
				 double dt,
				 double *phiOut,
				 double *wOut)
{
  int eN,i,id;
  double vol_eN,area_i,h_eN,d,volFact,areaFact;
  double nid,gphiNorm,gradHnorm,gradHdotGradNi,gradHdotGradPhi,tau,
    dhid,ai,gradNid;
  double gradPhi[3] = {0.,0.,0.};
  double K[4] = {0.,0.,0.,0.};
  double f[4] = {0.,0.,0.,0.};

  assert(nElementBoundaries_element == nSpace+1);
  assert(nDOF_element == nSpace+1);
  d = nSpace;
  /*mwf debug 
  printf("Entering advance LS SGS eN=%d neb=%d nde=%d Nsd=%d nq=%d nqb=%d\n",
	 nElements_global,nElementBoundaries_element,nDOF_element,
	 nSpace,nQuadraturePoints_element,nQuadraturePoints_elementBoundary);
  */
  /*convert metric tensors to volumes and areas for formulas*/
  volFact = 1.0; areaFact = 1.0;
  if (nSpace > 1)
    volFact = 0.5;
  if (nSpace > 2)
    { volFact = 1.0/6.0; areaFact = 0.5;}

  /*zero phiOut first*/
  for (eN = 0; eN < nElements_global; eN++)
    for (i=0; i < nDOF_element; i++)
      {
	phiOut[l2g[eN*nDOF_element+i]] = 0.0;
	wOut[l2g[eN*nDOF_element+i]] = 0.0;
      }
  for (eN = 0; eN < nElements_global; eN++)
    {
      h_eN = elementDiameters[eN];
      vol_eN = fabs(detJ[eN*nQuadraturePoints_element + 0])*volFact;/*assumed affine*/
      
      /*compute gradient of phiIn on element (constant)*/
      gradPhi[0] = 0.; gradPhi[1] = 0.; gradPhi[2] = 0.;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      /*nid is inner normal with area*/
	      nid = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradPhi[id] += phiIn[l2g[eN*nDOF_element+i]]*nid/(d*vol_eN);
	      /*mwf debug
	      printf("eN=%d ei=%i area=%g nid=%g gradPhi[%d]=%g\n",eN,i,area_i,nid,id,gradPhi[id]);
	      */
	    }/*end id*/
	  
	}/*i*/
      gphiNorm = 0.0;
      for (id = 0; id < nSpace; id++)
	gphiNorm += gradPhi[id]*gradPhi[id];
      gphiNorm = sqrt(gphiNorm);

      /*
	compute K[i] = \tau_i (\grad H_i \dot \grad N_i)(\grad H_i \dot \grad phi)
	relies on vertex quadrature used to eval dH and correspondence with loc dofs
      */
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/

	  gradHnorm = 0.0; gradHdotGradNi = 0.0; gradHdotGradPhi=0.0;
	  for (id = 0; id < nSpace; id++)
	    {
	      dhid = dH[eN*nQuadraturePoints_element*nSpace + i*nSpace + id];
	      nid  = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradNid = nid/(d*vol_eN);
	      
	      gradHnorm      += dhid*dhid;
	      gradHdotGradNi += dhid*gradNid;
	      gradHdotGradPhi+= dhid*gradPhi[id];

	    }
	  gradHnorm = sqrt(gradHnorm);
	  
	  /*tau formula from eqn (58) */
	  tau = sqrt(4.0/(dt*dt) + 4.0*gradHnorm*gradHnorm/(h_eN*h_eN) + 1.0e-12);
	  /*try 1 norm */
	  /*tau = 1./dt + 2.*gradHnorm/h_eN + 1.0e-12;*/
	  tau = 1.0/tau;

	  /*mwf debug
	  if (tau < 1.0e-4 || fabs(gradHdotGradNi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradNi=%g \n",eN,i,gradHdotGradNi);
	  if (fabs(gradHdotGradPhi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradPhi=%g \n",eN,i,gradHdotGradPhi);
	  */
	  /*mwf debug sharpens result for L_S_US_1d_circle 
	    tau = tau*0.01;
	  */
	  
	  K[i] = (1.0 + tau*gradHdotGradNi)*gradHdotGradPhi;
	  /*f = -r, vertex quad required*/
	  f[i] = -r[eN*nQuadraturePoints_element+i]*(1.+tau*gradHdotGradNi);
	  /*mwf debug
	  printf("SUPG lump eN=%d i=%d f[i]=%g zeroing\n",eN,i,f[i]);
	  */
	}/* i */
      /*final update*/
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  ai = vol_eN/(d+1.0);
	  phiOut[l2g[eN*nDOF_element+i]] += ai*(K[i]-f[i]); 
	  wOut[l2g[eN*nDOF_element+i]]   += ai;
	}
    }/*eN*/ 
  
}


void advanceStageRedistanceP1_C0_SUPG_lump(int nElements_global,
					   int nElementBoundaries_element,
					   int nDOF_element,
					   int nSpace,
					   int nQuadraturePoints_element,
					   int nQuadraturePoints_elementBoundary,
					   double* elementDiameters,
					   double* elementUnitOuterNormals,
					   double* detJ,
					   double* sqrtDetG,
					   double* elementQuadratureWeights,
					   int* l2g,
					   double * phiIn,
					   double * phiEvalIn,
					   double dt,
					   double eps,
					   double *phiOut,
					   double *wOut)
{
  int eN,i,id;
  double vol_eN,area_i,h_eN,d,volFact,areaFact;
  double nid,gphiNorm,gradHnorm,gradHdotGradNi,gradHdotGradPhi,tau,
    dhid,ai,gradNid,phiEval_i,He_i,S_i;
  double gradPhi[3] = {0.,0.,0.};
  double K[4] = {0.,0.,0.,0.};
  double f[4] = {0.,0.,0.,0.};

  assert(nElementBoundaries_element == nSpace+1);
  assert(nDOF_element == nSpace+1);
  d = nSpace;
  /*mwf debug
  printf("Entering advance LS SGS eN=%d neb=%d nde=%d Nsd=%d nq=%d nqb=%d\n",
	 nElements_global,nElementBoundaries_element,nDOF_element,
	 nSpace,nQuadraturePoints_element,nQuadraturePoints_elementBoundary);
  */
  /*convert metric tensors to volumes and areas for formulas*/
  volFact = 1.0; areaFact = 1.0;
  if (nSpace > 1)
    volFact = 0.5;
  if (nSpace > 2)
    { volFact = 1.0/6.0; areaFact = 0.5;}

  /*zero phiOut first*/
  for (eN = 0; eN < nElements_global; eN++)
    for (i=0; i < nDOF_element; i++)
      {
	phiOut[l2g[eN*nDOF_element+i]] = 0.0;
	wOut[l2g[eN*nDOF_element+i]] = 0.0;
      }
  for (eN = 0; eN < nElements_global; eN++)
    {
      h_eN = elementDiameters[eN];
      vol_eN = fabs(detJ[eN*nQuadraturePoints_element + 0])*volFact;/*assumed affine*/

      /*compute gradient of phiIn on element (constant)*/
      gradPhi[0] = 0.; gradPhi[1] = 0.; gradPhi[2] = 0.;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      /*nid is inner normal with area*/
	      nid = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradPhi[id] += phiIn[l2g[eN*nDOF_element+i]]*nid/(d*vol_eN);
	      /*mwf debug
	      printf("eN=%d ei=%i area=%g nid=%g gradPhi[%d]=%g\n",eN,i,area_i,nid,id,gradPhi[id]);
	      */
	    }/*end id*/
	  
	}/*i*/
      gphiNorm = 0.0;
      for (id = 0; id < nSpace; id++)
	gphiNorm += gradPhi[id]*gradPhi[id];
      gphiNorm = sqrt(gphiNorm);

      /*
	compute K[i] = \tau_i (\grad H_i \dot \grad N_i)(\grad H_i \dot \grad phi)
	relies on vertex quadrature used to eval dH and correspondence with loc dofs
      */

      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/

	  phiEval_i= phiEvalIn[l2g[eN*nDOF_element+i]];
	  He_i = 0.0;
	  if (fabs(phiEval_i) <= eps)
	    {
	      He_i = 0.5*(1.0+phiEval_i/eps + 1./M_PI*sin(M_PI*phiEval_i/eps));
	    }
	  else if (phiEval_i > eps)
	    He_i = 1.0;
	  S_i = 2.0*(He_i-0.5);

	  gradHnorm = 0.0; gradHdotGradNi = 0.0; gradHdotGradPhi=0.0;
	  for (id = 0; id < nSpace; id++)
	    {
	      dhid = S_i*gradPhi[id]/(gphiNorm+1.0e-10);
	      nid  = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradNid = nid/(d*vol_eN);
	      
	      gradHnorm      += dhid*dhid;
	      gradHdotGradNi += dhid*gradNid;
	      gradHdotGradPhi+= dhid*gradPhi[id];

	    }
	  gradHnorm = sqrt(gradHnorm);
	  
	  /*tau formula from eqn (58) */
	  tau = sqrt(4.0/(dt*dt) + 4.0*gradHnorm*gradHnorm/(h_eN*h_eN) + 1.0e-12);
	  /*try 1 norm */
	  /*tau = 1./dt + 2.*gradHnorm/h_eN + 1.0e-12;*/
	  tau = 1.0/tau;

	  /*mwf debug
	  if (tau < 1.0e-4 || fabs(gradHdotGradNi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradNi=%g \n",eN,i,gradHdotGradNi);
	  if (fabs(gradHdotGradPhi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradPhi=%g \n",eN,i,gradHdotGradPhi);
	  */
	  /*mwf debug sharpens result for L_S_US_1d_circle 
	    tau = tau*0.01;
	  */
	  
	  K[i] = (1.0 + tau*gradHdotGradNi)*gradHdotGradPhi;
	  /*f = -r, vertex quad required*/
	  f[i] = S_i*(1.+tau*gradHdotGradNi);
	  /*mwf debug
	  printf("test red eN=%d i=%d gradHdotGradNi=%g gradHdotGradPhi=%g ",eN,i
		 ,gradHdotGradNi,gradHdotGradPhi);
	  printf("tau=%g phiEval_i=%g Ki=%g fi=%g \n",tau,phiEval_i,K[i],f[i]);
	  */
	}/* i */
      /*final update*/
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  ai = vol_eN/(d+1.0);
	  phiOut[l2g[eN*nDOF_element+i]] += ai*(K[i]-f[i]); 
	  wOut[l2g[eN*nDOF_element+i]]   += ai;
	}
    }/*eN*/ 
}

void advanceStageRedistanceWeakDirP1_C0_SUPG_lump(int nElements_global,
						  int nElementBoundaries_element,
						  int nDOF_element,
						  int nSpace,
						  int nQuadraturePoints_element,
						  int nQuadraturePoints_elementBoundary,
						  int nExteriorElementBoundaries_global,
						  double* elementDiameters,
						  double* elementUnitOuterNormals,
						  double* detJ,
						  double* sqrtDetG,
						  double* elementQuadratureWeights,
						  int* l2g,
						  int* exteriorElementBoundariesArray,
						  int* elementBoundaryElementsArray,
						  double * phiIn,
						  double * phiEvalIn,
						  double dt,
						  double eps,
						  int* weakDirichletFlag,
						  double *phiOut,
						  double *wOut)
{
  int eN,i,id,ebN,extbN;
  double vol_eN,area_i,h_eN,d,volFact,areaFact;
  double nid,gphiNorm,gradHnorm,gradHdotGradNi,gradHdotGradPhi,tau,
    dhid,ai,gradNid,phiEval_i,He_i,S_i;
  double gradPhi[3] = {0.,0.,0.};
  double K[4] = {0.,0.,0.,0.};
  double f[4] = {0.,0.,0.,0.};
  double sign_eN=0.0; int lsCross_eN = 0;

  assert(nElementBoundaries_element == nSpace+1);
  assert(nDOF_element == nSpace+1);
  d = nSpace;
  /*mwf debug
  printf("Entering advance LS SGS eN=%d neb=%d nde=%d Nsd=%d nq=%d nqb=%d\n",
	 nElements_global,nElementBoundaries_element,nDOF_element,
	 nSpace,nQuadraturePoints_element,nQuadraturePoints_elementBoundary);
  */
  /*convert metric tensors to volumes and areas for formulas*/
  volFact = 1.0; areaFact = 1.0;
  if (nSpace > 1)
    volFact = 0.5;
  if (nSpace > 2)
    { volFact = 1.0/6.0; areaFact = 0.5;}

  /*zero phiOut first*/
  for (eN = 0; eN < nElements_global; eN++)
    for (i=0; i < nDOF_element; i++)
      {
	phiOut[l2g[eN*nDOF_element+i]] = 0.0;
	wOut[l2g[eN*nDOF_element+i]] = 0.0;
	weakDirichletFlag[l2g[eN*nDOF_element+i]] = 0;
      }
  for (eN = 0; eN < nElements_global; eN++)
    {
      h_eN = elementDiameters[eN];
      vol_eN = fabs(detJ[eN*nQuadraturePoints_element + 0])*volFact;/*assumed affine*/
      lsCross_eN = 0;
      sign_eN    = 0.0;

      /*compute gradient of phiIn on element (constant)*/
      gradPhi[0] = 0.; gradPhi[1] = 0.; gradPhi[2] = 0.;
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/
	  for (id = 0; id < nSpace; id++)
	    {
	      /*nid is inner normal with area*/
	      nid = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradPhi[id] += phiIn[l2g[eN*nDOF_element+i]]*nid/(d*vol_eN);
	      /*mwf debug
	      printf("eN=%d ei=%i area=%g nid=%g gradPhi[%d]=%g\n",eN,i,area_i,nid,id,gradPhi[id]);
	      */
	    }/*end id*/
	  
	}/*i*/
      gphiNorm = 0.0;
      for (id = 0; id < nSpace; id++)
	gphiNorm += gradPhi[id]*gradPhi[id];
      gphiNorm = sqrt(gphiNorm);

      /*
	compute K[i] = \tau_i (\grad H_i \dot \grad N_i)(\grad H_i \dot \grad phi)
	relies on vertex quadrature used to eval dH and correspondence with loc dofs
      */

      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  area_i = fabs(sqrtDetG[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary+
				 i*nQuadraturePoints_elementBoundary + 0])*areaFact; /*affine*/

	  phiEval_i= phiEvalIn[l2g[eN*nDOF_element+i]];
	  if (phiEval_i > 0.0 && sign_eN < 0.0)
	    lsCross_eN = 1;
	  else if (phiEval_i < 0.0 && sign_eN > 0.0)
	    lsCross_eN = 1;
	  else if (phiEval_i > 0.0 && sign_eN == 0.0)
	    sign_eN = 1.0;
	  else if (phiEval_i < 0.0 && sign_eN == 0.0)
	    sign_eN = -1.0;
	  else if (phiEval_i == 0.0)
	    {
	      lsCross_eN = 1; sign_eN = 1.0;
	    }
	  He_i = 0.0;
	  if (fabs(phiEval_i) <= eps)
	    {
	      He_i = 0.5*(1.0+phiEval_i/eps + 1./M_PI*sin(M_PI*phiEval_i/eps));
	    }
	  else if (phiEval_i > eps)
	    He_i = 1.0;
	  S_i = 2.0*(He_i-0.5);

	  gradHnorm = 0.0; gradHdotGradNi = 0.0; gradHdotGradPhi=0.0;
	  for (id = 0; id < nSpace; id++)
	    {
	      dhid = S_i*gradPhi[id]/(gphiNorm+1.0e-10);
	      nid  = -area_i
		*
		elementUnitOuterNormals[eN*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
					i*nQuadraturePoints_elementBoundary*nSpace+
					0*nSpace + id]; /*affine*/
	      /* \grad N_i = -|e_i|\vec n_i/d/|E| */
	      gradNid = nid/(d*vol_eN);
	      
	      gradHnorm      += dhid*dhid;
	      gradHdotGradNi += dhid*gradNid;
	      gradHdotGradPhi+= dhid*gradPhi[id];

	    }
	  gradHnorm = sqrt(gradHnorm);
	  
	  /*tau formula from eqn (58) */
	  tau = sqrt(4.0/(dt*dt) + 4.0*gradHnorm*gradHnorm/(h_eN*h_eN) + 1.0e-12);
	  /*try 1 norm */
	  /*tau = 1./dt + 2.*gradHnorm/h_eN + 1.0e-12;*/
	  tau = 1.0/tau;

	  /*mwf debug
	  if (tau < 1.0e-4 || fabs(gradHdotGradNi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradNi=%g \n",eN,i,gradHdotGradNi);
	  if (fabs(gradHdotGradPhi) < 1.0e-4)
	    printf("testImplSGS eN=%d i=%d gradHdotGradPhi=%g \n",eN,i,gradHdotGradPhi);
	  */
	  /*mwf debug sharpens result for L_S_US_1d_circle 
	    tau = tau*0.01;
	  */
	  
	  K[i] = (1.0 + tau*gradHdotGradNi)*gradHdotGradPhi;
	  /*f = -r, vertex quad required*/
	  f[i] = S_i*(1.+tau*gradHdotGradNi);
	  /*mwf debug
	  printf("test red eN=%d i=%d gradHdotGradNi=%g gradHdotGradPhi=%g ",eN,i
		 ,gradHdotGradNi,gradHdotGradPhi);
	  printf("tau=%g phiEval_i=%g Ki=%g fi=%g \n",tau,phiEval_i,K[i],f[i]);
	  */
	}/* i */
      /*final update*/
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  /*mwf force weak boundaries*/
	  if (lsCross_eN > 0)
	    weakDirichletFlag[l2g[eN*nDOF_element+i]] = 1;
	  ai = vol_eN/(d+1.0);
	  phiOut[l2g[eN*nDOF_element+i]] += ai*(K[i]-f[i]); 
	  wOut[l2g[eN*nDOF_element+i]]   += ai;
	}
    }/*eN*/ 
  /*fix up for freezing zero level sets*/
  /*first don't freeze on global boundary*/
  for (extbN= 0; extbN < nExteriorElementBoundaries_global; extbN++)
    {
      ebN = exteriorElementBoundariesArray[extbN];
      eN  = elementBoundaryElementsArray[ebN*2 + 0];
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  weakDirichletFlag[l2g[eN*nDOF_element+i]] = 0;
	  /*mwf debug
	  printf("supg red turning off weakDirFlag at boundary eN=%d i=%d \n",eN,i); 
	  */
	}
    } 
  for (eN = 0; eN < nElements_global; eN++)
    {
      for (i = 0; i < nElementBoundaries_element; i++)
	{
	  /* need to figure out how to avoid
 	     doing this on physical boundary*/
	  if (weakDirichletFlag[l2g[eN*nDOF_element+i]])
	    {
	      phiOut[l2g[eN*nDOF_element+i]] = 0.0;
	      /*mwf debug
	      printf("supg redist weakDirFlag eN=%d i=%d ig=%d \n",eN,i,l2g[eN*nDOF_element+i]);
	      */
	    }
	}
    }/*eN*/
}


/***********************************************************************
   Start DG limiting routines

 ***********************************************************************/
/* void computeElementAveragesP1Lagrange(int nElements_global, */
/* 				      int nDOF_element, */
/* 				      const int * l2g, */
/* 				      const double * Uin, */
/* 				      double * elementAverages) */
/* { */
/*   int eN,i,I; */
/*   const double denom = 1.0/nDOF_element; */
/*   for (eN = 0; eN < nElements_global; eN++) */
/*     { */
/*       elementAverages[eN] = 0.0; */
/*       for (i=0; i < nDOF_element; i++) */
/* 	{ */
/* 	  I = l2g[eN*nDOF_element+i]; */
/* 	  elementAverages[eN] += Uin[I]; */
/* 	} */
/*       elementAverages[eN] *= denom; */
/*     } */
/* } */
/* /\* */
/*   try simple sorting of A, but also keep track of sorted order */
/*   for other quantity based on A */
/* *\/ */
/* void shortSortDescending(double *A, int * order, int len) */
/* { */
/*   register int i,j,itmp; */
/*   register double tmp; */
/*   for (i=0; i < len; i++) */
/*     { */
/*       /\*order[i] = i;*\/ */
/*       for (j = i+1; j < len; j++) */
/* 	{ */
/* 	  if (A[j] > A[i]) */
/* 	    { */
/* 	      tmp = A[i]; A[i]= A[j];  A[j]= tmp; */
/* 	      itmp= order[i]; order[i] = order[j]; order[j] = itmp; */
/* 	    } */
/* 	} */
/*     } */
/* } */
/* void applyDurlofskyDGlimiterP1Lagrange3d(int killExtrema, */
/* 					 int allowMinWithUndershoot, */
/* 					 int nElements_global, */
/* 					 int nElementBoundaries_element, */
/* 					 int nNodes_element, */
/* 					 int nSpace, */
/* 					 int nDOF_element, */
/* 					 const int * elementNeighborsArray, */
/* 					 const int * elementBoundariesArray, */
/* 					 const int * elementNodesArray, */
/* 					 const double * nodeArray, */
/* 					 const double * elementBarycentersArray, */
/* 					 const double * elementBoundaryBarycentersArray, */
/* 					 const double * elementNeighborShapeGradients, */
/* 					 const int * l2g, */
/* 					 const double * grad_v0, */
/* 					 double * elementAverages, */
/* 					 int * tag, */
/* 					 double * Uin, */
/* 					 double * Uout) */
/* { */
/*   /\*try 3d version of  Durlofsky's original version with Phi limiter (not the min one)*\/ */
/*   int eN,nN_global,ebN,ebN_global,ebN_1,ebN_2,eN_ebN,eN_ebN_1,eN_ebN_2,onBoundary,isExtremum; */
/*   int dUdescendingOrder[5]; */
/*   double normDU[5],dU[5][3],max_ebN,min_ebN,uM; */
/*   /\*mwf for debugging*\/ */
/*   double normDUsave[5]; */
/*   int maxFound,minFound,islot,i,itmp,okGradient; */
/*   int nElementBoundaries_element2 = nElementBoundaries_element*nElementBoundaries_element; */
/*   register int DOF0,DOF1,DOF2,DOF3; */
/*   register double dUlim0,dUlim1,dUlim2; */
/*   register double uBar,uBar_ebN,uBar_ebN_1,uBar_ebN_2; */

/*   const double otol =  1.0e-5; */
/*   computeElementAveragesP1Lagrange(nElements_global,nDOF_element,l2g,Uin,elementAverages); */

/*   for (eN = 0; eN < nElements_global; eN++) */
/*     { */
/*       /\*current element*\/ */
/*       uBar = elementAverages[eN]; */
/*       DOF0 = l2g[eN*nDOF_element+0]; DOF1 = l2g[eN*nDOF_element+1];  */
/*       DOF2 = l2g[eN*nDOF_element+2]; DOF3 = l2g[eN*nDOF_element+3];  */
/*       onBoundary =  */
/* 	elementNeighborsArray[eN*nElementBoundaries_element + 0] < 0 || */
/* 	elementNeighborsArray[eN*nElementBoundaries_element + 1] < 0 || */
/* 	elementNeighborsArray[eN*nElementBoundaries_element + 2] < 0 || */
/* 	elementNeighborsArray[eN*nElementBoundaries_element + 3] < 0; */
/*       if (onBoundary > 0) */
/* 	{ */
/* 	  dUlim0 = 0.0; dUlim1 = 0.0; dUlim2 = 0.0; */
/* 	  tag[eN]  = 0; */
/* 	} */
/*       else */
/* 	{ */
/* 	  /\*check to see if this is a local extremum*\/ */
/* 	  maxFound = 0; minFound = 0; */
/* 	  for (ebN = 0; ebN < nElementBoundaries_element; ebN++) */
/* 	    { */
/* 	      eN_ebN        = elementNeighborsArray[eN*nElementBoundaries_element + ebN]; */
/* 	      uBar_ebN      = elementAverages[eN_ebN]; */
/* 	      maxFound = uBar_ebN > uBar ? 1 : maxFound; */
/* 	      minFound = uBar_ebN < uBar ? 1 : minFound; */
/* 	    } */
/* 	  isExtremum = (maxFound*minFound == 0); */
/* 	  if (isExtremum && killExtrema > 0) */
/* 	    { */
/* 	      /\*mwf debug */
/* 		printf("Durlofsky extremum found eN=%d uBar= %g uBar_ebN[%g,%g,%g]\n", */
/* 		eN,uBar,uBar_ebN[0],uBar_ebN[1],uBar_ebN[2]); */
/* 	      *\/ */
/* 	      dUlim0 = 0.0; dUlim1 = 0.0; dUlim2 = 0.0; */
/* 	      tag[eN]  = 0; */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      /\*by default take zero slope*\/ */
/* 	      dUlim0 = 0.0; dUlim1 = 0.0; dUlim2 = 0.0; */
/* 	      dU[0][0] =  */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 0*nSpace + 0]*Uin[DOF0]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 1*nSpace + 0]*Uin[DOF1]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 2*nSpace + 0]*Uin[DOF2]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 3*nSpace + 0]*Uin[DOF3]; */
/* 	      dU[0][1] =  */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 0*nSpace + 1]*Uin[DOF0]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 1*nSpace + 1]*Uin[DOF1]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 2*nSpace + 1]*Uin[DOF2]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 3*nSpace + 1]*Uin[DOF3]; */
/* 	      dU[0][2] =  */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 0*nSpace + 2]*Uin[DOF0]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 1*nSpace + 2]*Uin[DOF1]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 2*nSpace + 2]*Uin[DOF2]+ */
/* 		grad_v0[eN*1*nDOF_element*nSpace + 0*nDOF_element*nSpace + 3*nSpace + 2]*Uin[DOF3]; */
/* 	      normDU[0] = sqrt(dU[0][0]*dU[0][0] + dU[0][1]*dU[0][1] + dU[0][2]*dU[0][2]); */
/* 	      normDU[1] = 0.0; normDU[2] = 0.0; normDU[3] = 0.0; normDU[4] = 0.0; */
/* 	      /\*loop through neighbor simplexes and compute local interpolant gradients*\/ */
/* 	      for (ebN = 0; ebN < nElementBoundaries_element; ebN++) */
/* 		{ */
/* 		  ebN_1   = (ebN+1) % nElementBoundaries_element; */
/* 		  ebN_2   = (ebN+2) % nElementBoundaries_element; */
/* 		  /\*two neighbors for this simplex */
/* 		    local numbering is  */
/* 		    0 <--> this element */
/*                     1 <--> ebN neighbor */
/* 		    2 <--> ebN+1 neighbor */
/* 		    3 <--> ebN+2 neighbor */
/* 		   *\/ */
/* 		  eN_ebN    = elementNeighborsArray[eN*nElementBoundaries_element + ebN]; */
/* 		  eN_ebN_1  = elementNeighborsArray[eN*nElementBoundaries_element + ebN_1]; */
/* 		  eN_ebN_2  = elementNeighborsArray[eN*nElementBoundaries_element + ebN_2]; */
/* 		  uBar_ebN  = elementAverages[eN_ebN]; */
/* 		  uBar_ebN_1= elementAverages[eN_ebN_1]; */
/* 		  uBar_ebN_2= elementAverages[eN_ebN_2]; */

/* 		  /\*local number zero is always this element*\/ */
/* 		  dU[ebN+1][0] = 0.0; dU[ebN+1][1] = 0.0; dU[ebN+1][2] = 0.0; */
/* 		  dU[ebN+1][0]=  */
/* 		    uBar*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 						       ebN*nElementBoundaries_element*nSpace + */
/* 						       0*nSpace + 0] */
/* 		    + */
/* 		    uBar_ebN*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							   ebN*nElementBoundaries_element*nSpace + */
/* 							   1*nSpace + 0] */
/* 		    + */
/* 		    uBar_ebN_1*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							     ebN*nElementBoundaries_element*nSpace + */
/* 							     2*nSpace + 0] */
/* 		    + */
/* 		    uBar_ebN_2*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							     ebN*nElementBoundaries_element*nSpace + */
/* 							     3*nSpace + 0]; */
		  
/* 		  dU[ebN+1][1]=  */
/* 		    uBar*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 						       ebN*nElementBoundaries_element*nSpace + */
/* 						       0*nSpace + 1] */
/* 		    + */
/* 		    uBar_ebN*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							   ebN*nElementBoundaries_element*nSpace + */
/* 							   1*nSpace + 1] */
/* 		    + */
/* 		    uBar_ebN_1*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							     ebN*nElementBoundaries_element*nSpace + */
/* 							     2*nSpace + 1] */
/* 		    + */
/* 		    uBar_ebN_1*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							     ebN*nElementBoundaries_element*nSpace + */
/* 							     3*nSpace + 1]; */

/* 		  dU[ebN+1][2]=  */
/* 		    uBar*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 						       ebN*nElementBoundaries_element*nSpace + */
/* 						       0*nSpace + 2] */
/* 		    + */
/* 		    uBar_ebN*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							   ebN*nElementBoundaries_element*nSpace + */
/* 							   1*nSpace + 2] */
/* 		    + */
/* 		    uBar_ebN_1*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							     ebN*nElementBoundaries_element*nSpace + */
/* 							     2*nSpace + 2] */
/* 		    + */
/* 		    uBar_ebN_1*elementNeighborShapeGradients[eN*nElementBoundaries_element2*nSpace +  */
/* 							     ebN*nElementBoundaries_element*nSpace + */
/* 							     3*nSpace + 2]; */

/* 		  normDU[ebN+1] = sqrt(dU[ebN+1][0]*dU[ebN+1][0] + dU[ebN+1][1]*dU[ebN+1][1] + dU[ebN+1][2]*dU[ebN+1][2]); */
/* 		}/\*ebN*\/ */
/* 	      /\*now sort the gradients from largest to smallest*\/ */
/* 	      dUdescendingOrder[0] = 0; dUdescendingOrder[1] = 1; dUdescendingOrder[2]=2;  */
/* 	      dUdescendingOrder[3]=3;  dUdescendingOrder[4] = 4; */
/* 	      /\*save for debugging*\/ */
/* 	      normDUsave[0] = normDU[0]; normDUsave[1]=normDU[1]; normDUsave[2]=normDU[2]; normDUsave[3]=normDU[3]; */
/* 	      normDUsave[4] = normDU[4]; */
/* 	      shortSortDescending(normDU,dUdescendingOrder,5); */
/* 	      /\*mwf debug check ordering*\/ */
/* 	      for (i = 0; i < nElementBoundaries_element; i++) */
/* 		{ */
/* 		  if (normDU[i+1] > normDU[i]) */
/* 		    { */
/* 		      printf("\nPROBLEM Durlofsky 3d out of order eN=%d normDU[%d]=%g > normDU[%d] =%g \n", */
/* 			     eN,i+1,normDU[i+1], */
/* 			     i,normDU[i]); */
/* 		      printf("normDU=[%g,%g,%g,%g,%g] dUdescendingOrder=[%d,%d,%d,%d,%d]\n", */
/* 			     normDU[0],normDU[1],normDU[2],normDU[3],normDU[4], */
/* 			     dUdescendingOrder[0],dUdescendingOrder[1],dUdescendingOrder[2], */
/* 			     dUdescendingOrder[3],dUdescendingOrder[4]); */
/* 		      exit(1); */
/* 		    } */
/* 		  if (fabs(normDU[i]-normDUsave[dUdescendingOrder[i]]) > 1.0e-8) */
/* 		    { */
/* 		      printf("\nPROBLEM Durlofsky order wrong eN=%d normDU[%d]=%g normDUsave[%d] =%g \n", */
/* 			     eN,i,normDU[i],dUdescendingOrder[i],normDUsave[dUdescendingOrder[i]]); */
/* 		      printf("normDU=[%g,%g,%g,%g,%g] normDUsave=[%g,%g,%g,%g,%g] dUdescendingOrder=[%d,%d,%d,%d,%d]\n", */
/* 			     normDU[0],normDU[1],normDU[2],normDU[3],normDU[4], */
/* 			     normDUsave[0],normDUsave[1],normDUsave[2],normDUsave[3],normDUsave[4], */
/* 			     dUdescendingOrder[0],dUdescendingOrder[1],dUdescendingOrder[2],dUdescendingOrder[3], */
/* 			     dUdescendingOrder[4]); */
/* 		      exit(1); */

/* 		    } */
/* 		} */
	      
/* 	      /\*now start checking for overshoot, starting with largest du*\/ */
/* 	      okGradient = 0; islot = 0; */
/* 	      /\*use minimum slope if undershoot/overshoot detected for others*\/ */
/* 	      while (okGradient == 0 && islot < nElementBoundaries_element) */
/* 		{ */
/* 		  itmp = dUdescendingOrder[islot]; */
/* 		  /\*start with ok*\/ */
/* 		  okGradient = 1; */
/* 		  for (ebN = 0; ebN < nElementBoundaries_element; ebN++) */
/* 		    { */
/* 		      eN_ebN     = elementNeighborsArray[eN*nElementBoundaries_element + ebN]; */
/* 		      uBar_ebN   = elementAverages[eN_ebN]; */
/* 		      ebN_global = elementBoundariesArray[eN*nElementBoundaries_element + ebN]; */
/* 		      uM = uBar  */
/* 			+  */
/* 			dU[itmp][0]*(elementBoundaryBarycentersArray[ebN_global*3 + 0]- */
/* 				     elementBarycentersArray[eN*3 + 0]) */
/* 			+ */
/* 			dU[itmp][1]*(elementBoundaryBarycentersArray[ebN_global*3 + 1]- */
/* 				     elementBarycentersArray[eN*3 + 1]) */
/* 			+ */
/* 			dU[itmp][2]*(elementBoundaryBarycentersArray[ebN_global*3 + 2]- */
/* 				     elementBarycentersArray[eN*3 + 2]); */

/* 		      max_ebN = uBar > uBar_ebN ? uBar : uBar_ebN; */
/* 		      min_ebN = uBar < uBar_ebN ? uBar : uBar_ebN; */
		      
/* 		      okGradient = (okGradient > 0) && min_ebN <= uM -otol && uM <= max_ebN + otol; */
		      
/* 		    } */
/* 		  islot++;/\*undo this later if okGradient==1*\/ */
/* 		} */
/* 	      if (okGradient == 1) */
/* 		islot--; */
/* 	      if (okGradient || allowMinWithUndershoot > 0) */
/* 		{ */
/* 		  itmp = dUdescendingOrder[islot]; */
/* 		  dUlim0 = dU[itmp][0]; */
/* 		  dUlim1 = dU[itmp][1]; */
/* 		  dUlim2 = dU[itmp][2]; */
/* 		  tag[eN] = islot == 0 ? 1 : 0; */
/* 		  /\*mwf debug  */
/* 		    printf("Durlofsky eN=%d okGradient islot=%d dUlim[0]=%g dUlim[1]=%g dUlim[2]=%g\n", */
/* 		    eN,islot,dUlim0,dUlim1,dUlim2); */
/* 		  *\/ */
/* 		} */
	    
	      
/* 	    }/\*not an extremum*\/ */
/* 	}/\*in interior*\/ */
/*       /\*interpolate values to nodes assuming correspondence with local dofs*\/ */
/*       nN_global = elementNodesArray[eN*nNodes_element+0]; */
/*       Uout[DOF0] = uBar  */
/* 	+  */
/* 	dUlim0*(nodeArray[nN_global*3 + 0] - elementBarycentersArray[eN*3 + 0]) */
/* 	+ */
/* 	dUlim1*(nodeArray[nN_global*3 + 1] - elementBarycentersArray[eN*3 + 1]) */
/* 	+ */
/* 	dUlim2*(nodeArray[nN_global*3 + 2] - elementBarycentersArray[eN*3 + 2]); */
/*       nN_global = elementNodesArray[eN*nNodes_element+1]; */
/*       Uout[DOF1] = uBar  */
/* 	+  */
/* 	dUlim0*(nodeArray[nN_global*3 + 0] - elementBarycentersArray[eN*3 + 0]) */
/* 	+ */
/* 	dUlim1*(nodeArray[nN_global*3 + 1] - elementBarycentersArray[eN*3 + 1]) */
/* 	+ */
/* 	dUlim2*(nodeArray[nN_global*3 + 2] - elementBarycentersArray[eN*3 + 2]); */
/*       nN_global = elementNodesArray[eN*nNodes_element+2]; */
/*       Uout[DOF2] = uBar  */
/* 	+  */
/* 	dUlim0*(nodeArray[nN_global*3 + 0] - elementBarycentersArray[eN*3 + 0]) */
/* 	+ */
/* 	dUlim1*(nodeArray[nN_global*3 + 1] - elementBarycentersArray[eN*3 + 1]) */
/* 	+ */
/* 	dUlim2*(nodeArray[nN_global*3 + 2] - elementBarycentersArray[eN*3 + 2]); */
/*       nN_global = elementNodesArray[eN*nNodes_element+3]; */
/*       Uout[DOF3] = uBar  */
/* 	+  */
/* 	dUlim0*(nodeArray[nN_global*3 + 0] - elementBarycentersArray[eN*3 + 0]) */
/* 	+ */
/* 	dUlim1*(nodeArray[nN_global*3 + 1] - elementBarycentersArray[eN*3 + 1]) */
/* 	+ */
/* 	dUlim2*(nodeArray[nN_global*3 + 2] - elementBarycentersArray[eN*3 + 2]); */
      
/*     }/\*eN*\/ */
/* } */

/***********************************************************************
 try some different numerical fluxes
 ***********************************************************************/


/** @} */
