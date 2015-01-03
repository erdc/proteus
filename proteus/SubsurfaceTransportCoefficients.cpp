#include "SubsurfaceTransportCoefficients.h"
#include "pskRelations.h"
#include "densityRelations.h"

/*simple piecewise linear interpolation from a table
  assumes xv are increasing
 */
int findInterval(const double* vertices, int nv, double x, int* ival, double tol)
{
  int leftInt=0,rightInt=nv-2,failed=1,mid=0;
  assert(rightInt >= leftInt);
  /*take care of easy cases first*/
  if (fabs(x-vertices[leftInt]) < tol)
    {
      *ival=leftInt;
      failed=0;
      return failed;
    }
  if (x <= vertices[leftInt]-tol)
    {
      *ival=-1;
      failed=1;
      return failed;
    }
  if (fabs(x-vertices[rightInt+1]) < tol)
    {
      *ival=rightInt;
      failed=0;
      return failed;
    }
  if (x >= vertices[rightInt+1]+tol)
    {
      *ival = nv;
      failed=1;
      return failed;
    }
  /*otherwise, should have x in (left,right)*/
  while (leftInt <= rightInt)
    {
      mid = (int)(floor(0.5*(leftInt+rightInt)));
      if (vertices[mid] <= x && x < vertices[mid+1])/*x in interval mid*/
	{
	  *ival = mid;
	  failed = 0;
	  return failed;
	}
      else if (x < vertices[mid])/*x to the left of mid*/
	rightInt = mid-1;
      else if (x >= vertices[mid+1]) /*x to the right of mid*/
	leftInt = mid+1;
      else
	{
	  std::cout<<"findInterval shouldn't be here leftInt= "<<leftInt<<" rightInt= "<<rightInt<<std::endl;
	  assert(0);
	  failed = 1;
	  return failed;
	}
    }
  failed = 1;
  return failed;
} 
void piecewiseLinearTableLookup(double x,
				int nv,
				int* start,
				double* y,
				double* dy,
				const double* xv,
				const double* yv)
{
  int index=*start,findFailed=0;
  double val,tol=1.0e-8;
  findFailed = findInterval(xv,nv,x,&index,tol);
  if (findFailed && index == -1)
    {
      /*extrapolate off left, could use endpoint instead*/
      index=0;
    }
  else if (findFailed && index == nv)
    {
      /*extrapolate off right, could use endpoint instead*/
      index = nv-2;
    }
  else
    {
      assert(0 <= index && index < nv-1);
      assert(xv[index]-tol <= x && x<= xv[index+1]+tol);
    }
  assert(0 <= index && index < nv-1);
  val =  yv[index] +  (yv[index+1]-yv[index])/(xv[index+1]-xv[index])*(x-xv[index]);
  *y = val; 
  *dy = (yv[index+1]-yv[index])/(xv[index+1]-xv[index]); 
  *start = index;
}

#ifdef TRY_RE_NCP1
extern "C" int RE_NCP1_getResidual_VGM(double rho, //physical coefficients
				       double beta,
				       const double * gravity,
				       const double * alphaVG,
				       const double * nVG,
				       const double * thetaR,
				       const double * thetaSR,
				       const double * KWs,
				       const int * KWs_rowptr,
				       const int * KWs_colind,
				       //mesh info
				       int nSpace,
				       int nQuadraturePoints_element,
				       int nElements_global,
				       int nElementBoundaries_element,
				       int nElementBoundaries_global,
				       int nExteriorElementBoundaries_global,
				       int nQuadraturePoints_elementBoundary,
				       const double * nodeArray,
				       const int * elementNodesArray,
				       const double * elementBarycentersArray,
				       const int * elementMaterialTypes,
				       const int * elementNeighborsArray,
				       const int * elementBoundariesArray,
				       const double * elementBoundaryBarycentersArray,
				       const double * elementBoundaryOuterNormalsArray,
				       const int * exteriorElementBoundariesArray,
				       const int * elementBoundaryElementsArray,
				       const double * q_detJ,
				       //solution info
				       int nDOF_trial_element,
				       int nDOF_test_element,
				       const double * u_dof,
				       const double * u_l2g,
				       //time integration 
				       double alphaBDF,
				       const double * q_m_beta_bdf,
				       //storage to make upwinding simpler
				       double * q_kr,
				       //output for post-processing, time integration 
				       double * q_m,
				       double * q_u,
				       double * q_cfl,
				       double * q_elementResidual,
				       //boundary integral info, need to simplify
				       const int * isDOFBoundary_u,
				       const int * isFluxBoundary_u,
				       const double * ebqe_w_dS,
				       const double * ebqe_u,
				       const double * ebqe_flux
				       //residual info
				       int offset_u,
				       int stride_u,
				       double * globalResidal,
				       )
{
  assert (nSpace+1 == nQuadraturePoints_element);
  assert (nElementBoundaries_element == nSpace+1);
  assert (nDOF_trial_element == nSpace+1);
  //temporaries
  double krw,dkrw,thetaw,dthetaw,rhom;
  double KWs_ebN[9];
  const double nAvgWeight = 1.0/(nSpace+1.);
  double volFactor = 1.0;
  if (nSpace == 2)
    volFactor = 0.5;
  if (nSpace == 3)
    volFactor = 1.0/6.0;
  const int nnz = KWs_rowptr[nSpace];
  double w[4],grad_w[4][3],grad_u[3];

  //evaluate mass term and element psk relations in first loop
  //store kr for upwinding, then loop again for spatial integrals  
  for (int eN = 0; eN < nElements_global; eN++)
    {
      const double volume = volFactor * fabs(q_detJ[eN*nQuadraturePoints_element]);//affine
      const int matID = elementMaterialTypes[eN];
      const double weight = nAvgWeight*volume;

      for (int ebN=0; ebN < nElementBoundaries_element; ebN++)
	{
	  const int eN_ebN = eN*nDOF_test_element+ebN;
	  const double u_j = u_dof[u_l2g[eN_ebN]];
	  const double psic = -u_j;
	  pskEvaluate_vgm_re(psic,
			     alphaVG[matID],
			     nVG[matID],
			     thetaR[matID],
			     thetaSR[matID],
			     krw,thetaw,
			     dkrw,dthetaw);
	  
	  q_kr[eN_ebN] = krw;
	  q_u[eN_ebN]  = u_j;
	  //(diagonal) mass term
	  rhom = rho*exp(beta*u_j);
	  //save for time history
	  q_m[ebN_ebN] = rhom*thetaw;
	  const double m_t_ebN = alphaBDF*q_m + q_m_beta_bdf[eN_ebN];

	  q_elementResidual[eN_ebN] += weight*m_t_ebN;
	  globalResidual[offset_u+stride_u*u_l2g[eN_ebN]] += weight*m_t_ebN;

	}//ebN
    }//eN

  //stiffness term
  for (int eN = 0; eN < nElements_global; eN++)
    {
      evaluateTestFunctionsAndGradientsOnElement(nSpace,
						 eN,
						 &elementBarycentersArray[eN*3],
						 nodeArray,
						 elementBoundaryOuterNormalsArray,
						 w,
						 grad_w);
      //constant on an element
      for (int I = 0; I < nSpace; I++)
	{
	  grad_u[I] = 0.0;
	  for (int j = 0; j < nDOF_trial_element; j++)
	    {
	      grad_u[I] += u_dof[u_l2g[eN*nDOF_trial_element+j]]*grad_w[j][I];
	    }
	}
      double kr_eN = 0.0,u_eN=0.0;
      //averages for kr, upwinding
      for (int ebN=0; ebN < nElementBoundaries_element; ebN++)
	{
	  const int eN_ebN = eN*nDOF_test_element+ebN;
	  kr_eN += q_kr[eN_ebN]*nAvgWeight;
	  u_eN  += q_u[eN_ebN]*nAvgWeight;
	}
      //for upwinding
      const double phi_eN = u_eN;
      for (int I = 0; I < nSpace; I++)
	phi_eN -= gravity[I]*elementBarycentersArray[eN*3+I];
      //loop over each element, determine their contribution to each residual 
      //equation's stiffness term
      for (int ebN=0; ebN < nElementBoundaries_element; ebN++)
	{
	  const int eN_neighbor = elementNeighborsArray[eN*nElementBoundaries_element + ebN];
	  //by default current element is upwind
	  double kr_up = kr_eN; 
	  if (eN_neighbor >= 0)
	    {
	      const matID_neighbor = elementMaterialTypes[eN_neighbor];
	      //harmonic average for intrinsic perm
	      for (int ii=0; ii < nnz; ii++)
		KWs_ebN[ii] = 2.0*KWs[matID*nnz+ii]*KWs[matID_neighbor*nnz+ii]/(KWs[matID*nnz+ii] + KWs[matID_neighbor*nnz+ii] + 1.0e-20);
	      
	      //averages for kr, upwinding
	      double kr_neig = 0.0,u_neig=0.0;

	      for (int ebN_neig=0; ebN_neig < nElementBoundaries_element; ebN_neig++)
		{
		  const int eN_ebN_neig = eN_neighbor*nDOF_test_element+ebN_neig;
		  kr_neig += q_kr[eN_ebN_neig]*nAvgWeight;
		  u_neig  += q_u[eN_ebN_neig]*nAvgWeight;
		}
	    
	      const double phi_neig = u_neig;
	      for (int I = 0; I < nSpace; I++)
		phi_neig -= gravity[I]*elementBarycentersArray[eN_neighbor*3+I];
	      
	      if (phi_neig > phi_eN)
		{
		  kr_up = kr_neig;
		}
	    }//interior face
	  else //exterior face
	    {
	      for (int ii=0; ii < nnz; ii++)
		KWs_ebN[ii] = KWs[matID*nnz+ii];
	      //todo could upwind using Dirichlet bc
	    }
	  //loop through element residuals and accumulate contribution of this face to advection and stiffness integrals
	  for (int i =0; i < nDOF_test_element; i++)
	    {
	      const int eN_i = eN*nDOF_test_element + i;
	      for (int I = 0; I < nSpace; I++)
		{
		  for (int ii=KWs_rowptr[I]; ii < KWs_rowptr[I+1]; ii++)
		    {
		      q_elementResidual[eN_i] -= rho*rho*KWs_ebN[ii]*gravity[KWs_colind[ii]]*grad_w[i][I];
		      globalResidual[offset_u+stride_u*u_l2g[eN_i]] -= rho*rho*KWs_ebN[ii]*gravity[KWs_colind[ii]]*grad_w[i][I];
		      q_elementResidual[eN_i] += rho*KWs_ebN[ii]*grad_u[KWs_colind[ii]]*grad_w[i][I];
		      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += rho*KWs_ebN[ii]*grad_u[KWs_colind[ii]]*grad_w[i][I];

		    }//ii
		}//I
	    }//test functions
	}//element boundaries
    }//eN stiffness integrals

  //boundary conditions, Dirichlet done strongly
  //could assume ebN agrees with global dof, 
  //but need q_elementResidual for post-processing
  //could use ebN --> eN,ebN_local mapping to get which test function is
  //nonzero on the element
  //todo need to decide about how to get element boundary areas
  double elementResidual[4];
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      const int ebN = exteriorElementBoundariesArray[ebNE];
      const int eN  = elementBoundaryElementsArray[ebN*2+0];

      for (int kb = 0; kb < nQuadraturePoints_elementBoundary; kb++)
	{
	  const int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary + kb;
	  if (isFluxBoundary_u[ebNE_kb])
	    {
	      for (int i = 0; i < nDOF_test_element; i++)
		{
		  elementResidual[i] += ebqe_flux[ebNE_kb]*ebqe_w_dS[ebNE_kb];
		}
	    }
	}//kb
      for (int i = 0; i < nDOF_test_element; i++)
	{
	  q_elementResidual[eN*nDOF_test_element + i] += elementResidual[i];
	  globalResidual[offset_u+stride_u*u_l2g[eN*nDOF_test_element+i]] += elementResidual[i];
	}
    }//ebNE

  //apply dirichlet bc's strongly?
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      const int ebN = exteriorElementBoundariesArray[ebNE];
      const int eN  = elementBoundaryElementsArray[ebN*2+0];
      //assume bc's constant on a face
      if (isDOFBoundary_u[ebNE*nQuadraturePoints_elementBoundary])
	{
	  for (int j = 0; j < nElementBoundaries_element; j++)
	    {
	      if (elementBoundariesArray[eN*nElementBoundaries_element + j] == ebN)
		{
		  globalResidual[offset_u+stride_u*u_l2g[eN*nDOF_trial_element+j]] = u_dof[u_l2g[eN*nDOF_trial_element+j]] - ebqe_u[ebNE*nQuadraturePoints_elementBoundary];
		}
	    }
	}
    }
}
//end TRY_RE_NCP1
#endif 
//--------------------------------------------------
//calculate Rusanov flux at global boundaries using
//piecewise constant information
//--------------------------------------------------
int calculateRusanovFluxSaturationEquationIncomp_PWC(double safetyFactor,
						     int nSpace,
						      //physical information
						      //psk type
						      int pskModelFlag,
						      int nParams,
						      const int* rowptr,
						      const int* colind,
						      const int* materialTypes,
						      double muw,
						      double mun,
						      const double* omega,
						      const double* Kbar,
						      double b,
						      const double* rwork_psk,
						      const double* rwork_psk_tol,
						      const double* rwork_density_w,
						      const double* rwork_density_n,
						      const double* g,
						      //velocity at interface 
						      const double* ebq_global_qt,
						      //bounds on characteristic speed for elements
						      const double* q_lambda_bar,
						      //mesh information
						      int nElements_global,
						      int nElementBoundaries_element,
						      int nInteriorElementBoundaries_global,
						      int nExteriorElementBoundaries_global,
						      int nQuadraturePoints_element, 
						      int nQuadraturePoints_elementBoundary,
						      const int* interiorElementBoundariesArray,
						      const int* exteriorElementBoundariesArray,
						      const int* elementBoundaryElementsArray,
						      const int* elementBoundaryLocalElementBoundariesArray,
						      const double* n,
						      //solution information
						      const double * q_u,
						      //element quadrature
						      const double* dV,
						      //int nDOF_trial_element,
						      //const int* u_l2g,
						      //const double* u_dof,
						      //boundary information
						      const int * isDOFBoundary,
						      const double* bc_u,
						      //output
						      double* flux)
{
  //loop over interior faces, compute average solution value for right and left
  //assumes characteristic speed bound has been computed on each element at 
  //quadrature points 
  //assuming incompressible fractional flow formulation
  //apply rusanov flux formula to calculate numerical flux at face assuming oriented
  //in direction of 'left' aka 0 element neighbor
  //
  //loop over exterior faces, repeat averaging step and
  //apply rusanov flux at dof boundaries using bc condition assumed constant over face
  //as 'right' state

  int failed =0;
  //coefficients for 2-phase flow incompressible (slightly compressible should be same)
  PskRelation* psk;
  switch (pskModelFlag)
    {
    case 0:
      psk = new SimplePSK(rwork_psk);
      break;
    case 1:
      psk = new VGM(rwork_psk);
      break;
    case 2:
      psk = new VGB(rwork_psk);
      break;
    case 3:
      psk = new BCM(rwork_psk);
      break;
    case 4:
      psk = new BCB(rwork_psk);
      break;
    default:
      std::cerr<<"pskModelFlag= "<<pskModelFlag<<" not recognized"<<std::endl;
      assert (false);
      break;
    }//psk model selection
  psk->setTolerances(rwork_psk_tol);
  FractionalFlowVariables fracFlow(muw,mun);
  //normalized densities \rho_{\alpha} = \varrho_{\alpha}/\varrho_{\alpha,0}
  //  for spatial term, assuming slight compressiblity so assume \rho_{\alpha} = 1
  const double rwork_density_w_x[1] = {1.0}; const double rwork_density_n_x[1] = {1.0}; 
  ConstantDensity density_w_x(rwork_density_w_x),density_n_x(rwork_density_n_x);
  const int nnz = rowptr[nSpace];
  double f_left[3] = {0.,0.,0.}; double f_right[3] = {0.,0.,0.};

  for (int ebNI = 0; ebNI < nInteriorElementBoundaries_global; ebNI++)
    {
      const int ebN = interiorElementBoundariesArray[ebNI];
      const int eN_left  =elementBoundaryElementsArray[ebN*2+0];
      const int eN_right =elementBoundaryElementsArray[ebN*2+1];
      const int ebN_left =elementBoundaryLocalElementBoundariesArray[ebN*2+0];
      const int ebN_right=elementBoundaryLocalElementBoundariesArray[ebN*2+1];

      //compute:
      //  left and right states
      double u_left=0.0,u_right=0.0,vol_left=0.0,vol_right=0.0;
      //  max speeds on element
      double maxSpeed_element=0.,maxSpeed_left=0.0,maxSpeed_right=0.0,tmp_left=0.0,tmp_right=0.0;

      for (int k = 0; k < nQuadraturePoints_element; k++)
	{
	  u_left   += q_u[eN_left*nQuadraturePoints_element+k]*dV[eN_left*nQuadraturePoints_element+k];
	  vol_left += dV[eN_left*nQuadraturePoints_element+k];
	  u_right  += q_u[eN_right*nQuadraturePoints_element+k]*dV[eN_right*nQuadraturePoints_element+k];
	  vol_right+= dV[eN_right*nQuadraturePoints_element+k];
	  tmp_left  = q_lambda_bar[eN_left*nQuadraturePoints_element+k];
	  tmp_right = q_lambda_bar[eN_right*nQuadraturePoints_element+k];
	  maxSpeed_left = fabs(tmp_left)  > maxSpeed_left ? fabs(tmp_left) : maxSpeed_left;
	  maxSpeed_right= fabs(tmp_right) > maxSpeed_right ? fabs(tmp_right) : maxSpeed_right;	  
	}
      u_left /= vol_left;
      u_right/= vol_right;
      maxSpeed_element = maxSpeed_right > maxSpeed_left ? maxSpeed_right : maxSpeed_left;

      //compute fluxes at left and right state assuming two-phase saturation equation
      //since just using averages compute values and save them
      //left
      const int matID_left = materialTypes[eN_left];
      psk->setParams(&rwork_psk[matID_left*nParams]);
      psk->calc(u_left);
      fracFlow.calc(*psk,density_w_x,density_n_x);
      const double fw_left=fracFlow.fw; const double fn_left=fracFlow.fn; const double law_left=fracFlow.lambdaw;
      const double rho_n_left=density_n_x.rho; const double rho_w_left = density_w_x.rho;
      //right
      const int matID_right = materialTypes[eN_right];
      psk->setParams(&rwork_psk[matID_right*nParams]);
      psk->calc(u_right);
      fracFlow.calc(*psk,density_w_x,density_n_x);
      const double fw_right=fracFlow.fw; const double fn_right=fracFlow.fn; const double law_right=fracFlow.lambdaw;
      const double rho_n_right=density_n_x.rho; const double rho_w_right = density_w_x.rho;

      //stuck with quadrature point loop for now
      for (int k=0; k < nQuadraturePoints_elementBoundary; k++)
	{
	  //redundant except for qt
	  for (int I=0; I < nSpace; I++)
	    {
	      f_left[I] = ebq_global_qt[ebN*nQuadraturePoints_elementBoundary*nSpace + k*nSpace + I]*fw_left;
	      //todo think through sign
	      f_right[I]= ebq_global_qt[ebN*nQuadraturePoints_elementBoundary*nSpace + k*nSpace + I]*fw_right;
	      for (int m=rowptr[I]; m < rowptr[I+1]; m++)
		{
		  const int J = colind[m];
		  f_left[I] -= Kbar[matID_left*nnz+m]*law_left*fn_left*(b*rho_n_left-rho_w_left)*g[J]; 
		  f_right[I]-= Kbar[matID_right*nnz+m]*law_right*fn_right*(b*rho_n_right-rho_w_right)*g[J]; 
		}
	    }

	  double left_flux=0.0,right_flux=0.0;
	  for (int I=0; I < nSpace; I++)
	    {
	      left_flux += n[eN_left*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
			     ebN_left*nQuadraturePoints_elementBoundary*nSpace+k*nSpace+I]
		*
		f_left[I];
	      right_flux += n[eN_left*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
			      ebN_left*nQuadraturePoints_elementBoundary*nSpace+k*nSpace+I]
		*
		f_right[I];
	    }

	  const double maxSpeed =maxSpeed_element*safetyFactor;
	  flux[ebN*nQuadraturePoints_elementBoundary+k] = 0.5*(left_flux+right_flux) -0.5*maxSpeed*(u_right-u_left);
	}//k

    }//ebNI

  //exterior element boundaries
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      const int ebN = exteriorElementBoundariesArray[ebNE];
      const int eN_left  =elementBoundaryElementsArray[ebN*2+0];
      const int ebN_left =elementBoundaryLocalElementBoundariesArray[ebN*2+0];

      //compute:
      //  left and right states
      double u_left=0.0,u_right=0.0,vol_left=0.0;
      //  max speeds on element
      double maxSpeed_element=0.,maxSpeed_left=0.0,maxSpeed_right=0.0,tmp_left=0.0,tmp_right=0.0;

      for (int k = 0; k < nQuadraturePoints_element; k++)
	{
	  u_left   += q_u[eN_left*nQuadraturePoints_element+k]*dV[eN_left*nQuadraturePoints_element+k];
	  vol_left += dV[eN_left*nQuadraturePoints_element+k];
	  tmp_left  = q_lambda_bar[eN_left*nQuadraturePoints_element+k];
	  maxSpeed_left = fabs(tmp_left)  > maxSpeed_left ? fabs(tmp_left) : maxSpeed_left;
	}
      u_left /= vol_left;

      //todo compute bounds for speed using dirichlet information
      maxSpeed_element = maxSpeed_left;

      //compute fluxes at left and right state assuming two-phase saturation equation
      //since just using averages compute values and save them
      //left
      const int matID_left = materialTypes[eN_left];
      psk->setParams(&rwork_psk[matID_left*nParams]);
      psk->calc(u_left);
      fracFlow.calc(*psk,density_w_x,density_n_x);
      const double fw_left=fracFlow.fw; const double fn_left=fracFlow.fn; const double law_left=fracFlow.lambdaw;
      const double rho_n_left=density_n_x.rho; const double rho_w_left = density_w_x.rho;


      //right
      const int matID_right = materialTypes[eN_left];
      psk->setParams(&rwork_psk[matID_right*nParams]);

      //stuck with quadrature point loop for now
      for (int k=0; k < nQuadraturePoints_elementBoundary; k++)
	{
	  //redundant except for qt
	  for (int I=0; I < nSpace; I++)
	    {
	      f_left[I] = ebq_global_qt[ebN*nQuadraturePoints_elementBoundary*nSpace + k*nSpace + I]*fw_left;
	      for (int m=rowptr[I]; m < rowptr[I+1]; m++)
		{
		  const int J = colind[m];
		  f_left[I] -= Kbar[matID_left*nnz+m]*law_left*fn_left*(b*rho_n_left-rho_w_left)*g[J]; 
		}
	    }
	  if (isDOFBoundary[ebNE*nQuadraturePoints_elementBoundary+k])
	    {
	      u_right = bc_u[ebNE*nQuadraturePoints_elementBoundary+k];
	      psk->calc(u_right);
	      fracFlow.calc(*psk,density_w_x,density_n_x);
	      const double fw_right=fracFlow.fw; const double fn_right=fracFlow.fn; const double law_right=fracFlow.lambdaw;
	      const double rho_n_right=density_n_x.rho; const double rho_w_right = density_w_x.rho;
	      
	      for (int I=0; I < nSpace; I++)
		{
		  //todo think through sign
		  f_right[I]= ebq_global_qt[ebN*nQuadraturePoints_elementBoundary*nSpace + k*nSpace + I]*fw_right;
		  for (int m=rowptr[I]; m < rowptr[I+1]; m++)
		    {
		      const int J = colind[m];
		      f_right[I]-= Kbar[matID_right*nnz+m]*law_right*fn_right*(b*rho_n_right-rho_w_right)*g[J]; 
		    }
		}
	    }
	  else
	    {
	      u_right=u_left;
	      for (int I=0; I < nSpace; I++)
		f_right[I] = f_left[I];
	    }

	  double left_flux=0.0,right_flux=0.0;
	  for (int I=0; I < nSpace; I++)
	    {
	      left_flux += n[eN_left*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
			     ebN_left*nQuadraturePoints_elementBoundary*nSpace+k*nSpace+I]
		*
		f_left[I];
	      right_flux += n[eN_left*nElementBoundaries_element*nQuadraturePoints_elementBoundary*nSpace+
			      ebN_left*nQuadraturePoints_elementBoundary*nSpace+k*nSpace+I]
		*
		f_right[I];
	    }

	  const double maxSpeed =maxSpeed_element*safetyFactor;
	  flux[ebN*nQuadraturePoints_elementBoundary+k] = 0.5*(left_flux+right_flux) -0.5*maxSpeed*(u_right-u_left);
	}//k

    }//ebNE

  //clean up
  delete psk;

  return failed;
}

