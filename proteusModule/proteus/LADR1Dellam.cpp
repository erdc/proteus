#include "LADR1Dellam.h"
#include "tracking.h"

/***********************************************************************

  TODO:
       old mass integral
          break assumption about x_track having same layout as q_x,q_u?
       source term integral
       outflow boundary condition
          see if need to debug/improve approx. for large Cr
       full tensor diffusion coefficient
       clean up const declarations 
       add linear varying velocity approximation for tracking
 **********************************************************************/


extern "C" void calculateResidual_LADR1Dellam(double theta,    //physical parameters
					      double aL,
					      double Dm,
					      double dtnp1,                 //full time step
					      int nElements_global,//mesh representation
					      int nNodes_global,
					      int nNodes_element,
					      int nElementBoundaries_element,
					      int nQuadraturePoints_element,
					      const double * nodeArray,
					      const int * elementNodesArray,
					      const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
					      const double* dV,             //integration weights 
					      const double* x_track,        //location of forward tracked integration points, assumed for now 
                                                                      //nElements_global x nQuadraturePoints_element
					      const double* t_track,        //time forward tracked points stopped (t^{n+1} or earlier if exited domain)
					      const int* element_track,     //element each forward tracked point ended up in
					      const int* flag_track,        //id for each point, 0 -- interior, 1 exited domain, -1 didn't track for some reason
					      double* x_dt,                 //spatially varying time steps
					      int* u_l2g, 
					      double* elementDiameter,
					      double* u_dof,
					      double* u_trial, 
					      double* u_grad_trial, 
					      double* u_test_dV, 
					      double* u_grad_test_dV, 
					      double* velocity,
					      double* q_u,
					      double* q_m,
					      double* q_m_last,
					      double* q_cfl,
					      double* q_elementResidual_u, 
					      int* sdInfo_u_rowptr, int* sdInfo_u_colind,
					      int offset_u, int stride_u, 
					      double* globalResidual,
					      int nExteriorElementBoundaries_global,
					      int nQuadraturePoints_elementBoundary,
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
					      double* ebqe_outflow_flux_last,
					      double* u_test_dS_ext,
					      double* ebqe_u,
					      double* ebqe_outflow_flux)//save outflow flux for trapezoidal rule approximation
{
  using namespace LADR1Dellam;
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
	    m=0.0,dm=0.0,m_old=0,
	    f[nSpace],df[nSpace],a[nSpace],//todo make full tensor
	    dm_t=0.0;
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
	  //todo make a full tensor add aT
          evaluateCoefficients_c(theta,
				 &velocity[eN_k_nSpace],//make sure this has q in it
				 aL,
				 Dm,
				 u,
				 m,
				 dm,
				 f,
				 df,
				 a);
	  dm_t  = dm/x_dt[eN_k];
	  //
          //calculte cfl still for now, could calculate Peclet number too
	  //
	  calculateCourantNumber_c(elementDiameter[eN],dm,df,q_cfl[eN_k]);

          // 
          //update element residual 
          // 
          for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_k_i=eN_k*nDOF_test_element+i,
		eN_k_i_nSpace = eN_k_i*nSpace;
	      //mwf debug
	      //std::cout<<"ellam resid eN= "<<eN<<" k= "<<k<<" m= "<<m<<" dt= "<<x_dt[eN_k]<<std::endl;
	      elementResidual_u[i] += Mass_weak_c(m,u_test_dV[eN_k_i]) + 
  		Diffusion_weak_c(x_dt[eN_k],sdInfo_u_rowptr,sdInfo_u_colind,a,grad_u,&u_grad_test_dV[eN_k_i_nSpace]);
            }//i
	  //
	  //save solution for other models and time integration 
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
          globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
        }//i
    }//elements


  //
  //calculate old mass contribution
  //
  /***********************************************************************
    for current integration point x_k
       1. get mass at old time level, m^{n}_k
       2. get integration weight W_k
       3. get forward tracked integration point x^{n+1}_k
       4. evaluate test functions with non-zero support at x^{n+1}_k
           a. get element eN^{n+1}_k from element_track[k]
           b. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       5. accumulate m^{n+1}_k*W_k*w_i(x^{n+1}_{k}) to global residual I(i) 
  **********************************************************************/
  //todo switch to a loop over global integration points?
  for(int eN=0;eN<nElements_global;eN++)
    {
      //declare local storage for element residual and initialize
      register double w[nDOF_test_element];
      //loop over quadrature points, may end up on elements with different support
      for  (int k=0;k<nQuadraturePoints_element;k++)
        {
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      w[i] = 0.0;
	    }
	  //compute indeces and declare local storage
	  register int eN_k = eN*nQuadraturePoints_element+k;
	  //todo decide if ignoring outflow mass is a better idea
	  if (flag_track[eN_k] >= 0)
	    {
	      //m_k^{n+1},W_k
	      register double m_old = q_m_last[eN_k],
		weight = dV[eN_k];
	      register int eN_track = element_track[eN_k];
	      //mwf debug
	      //std::cout<<"ellam tracked point ["<<eN<<","<<k<<"] --> "<<x_track[eN_k*3]
	      //	       <<" eN_track= "<<eN_track<<" x0= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3]
	      //       <<" x1= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3]<<std::endl;
	      evaluateTestFunctionsOnElement(&x_track[eN_k*3],
					     &nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3],
					     &nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3],
					     w);
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int eN_track_i=eN_track*nDOF_test_element+i;
		  q_elementResidual_u[eN_track_i] -= m_old*w[i]*weight;
		  globalResidual[offset_u+stride_u*u_l2g[eN_track_i]] -= m_old*w[i]*weight; 
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element
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
	    a_ext[nSpace],//todo make full tensor
	    flux_ext=0.0,
	    bc_u_ext=0.0,
	    bc_grad_u_ext[nSpace],
	    bc_m_ext=0.0,
	    bc_dm_ext=0.0,
	    bc_f_ext[nSpace],
	    bc_df_ext[nSpace],
	    bc_a_ext[nSpace];//todo make full tensor
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
	  evaluateCoefficients_c(theta,
				 &ebqe_velocity_ext[ebNE_kb_nSpace],
				 aL,
				 Dm,
				 u_ext,
				 m_ext,
				 dm_ext,
				 f_ext,
				 df_ext,
				 a_ext);
	  evaluateCoefficients_c(theta,
				 &ebqe_velocity_ext[ebNE_kb_nSpace],
				 aL,
				 Dm,
				 bc_u_ext,
				 bc_m_ext,
				 bc_dm_ext,
				 bc_f_ext,
				 bc_df_ext,
				 bc_a_ext);    
	  //save for other models?
	  ebqe_u[ebNE_kb] = u_ext;
	  // 
	  //calculate outflow numerical fluxes, 
	  //for now applies zero diffusive flux if v.n >= 0 
	  // 
	  exteriorOutflowFlux_c(&ebqe_n_ext[ebNE_kb_nSpace],
				u_ext,
				&ebqe_velocity_ext[ebNE_kb_nSpace],
				flux_ext);
	  //save for next time step
	  ebqe_outflow_flux[ebNE_kb] = flux_ext;
	  //add time scaling
	  flux_ext = dtnp1*0.5*(flux_ext + ebqe_outflow_flux_last[ebNE_kb]);
	  //mwf debug
	  //std::cout<<"ebNE= "<<ebNE<<" u_ext= "<<u_ext<<" ebqe_bc_u_ext= "<<ebqe_bc_u_ext[ebNE_kb] <<" bc_u_ext= "<<bc_u_ext
	  //	   <<" isDOFBoundary_u= "<<isDOFBoundary_u[ebNE_kb]<<" v[0]= "<<ebqe_velocity_ext[ebNE_kb_nSpace]
	  //	   <<" isFluxBoundary_u= "<<isFluxBoundary_u[ebNE_kb]<<" flux_ext= "<<flux_ext
	  //	   <<" flux_ext_last= "<<ebqe_outflow_flux_last[ebNE_kb]<<std::endl;

	  //
	  //update residuals
	  //
	  //
	  //todo make sure scale flux by correct (full dt) in boundary approximation
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;

	      elementResidual_u[i] += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
	      //globalConservation += ExteriorElementBoundaryFlux_c(flux_ext,u_test_dS_ext[ebNE_kb_i]);
	    }//i
	}//kb
      //
      //update the element and global residual storage
      //
      for (int i=0;i<nDOF_test_element;i++)
	{
	  int eN_i = eN*nDOF_test_element+i;

	  q_elementResidual_u[eN_i] += elementResidual_u[i];
	  globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
	}//i
    }//ebNE
  //std::cout<<"Ellam global conservation============================================================="<<globalConservation<<std::endl;
}


extern "C" void calculateJacobian_LADR1Dellam (double theta,    //physical parameters
					       double aL,
					       double Dm,
					       double dtnp1,
					       int nElements_global,
					       int nQuadraturePoints_element,
					       double* x_dt,
					       int* u_l2g,
					       double* elementDiameter,
					       double* u_dof, 
					       double* u_trial, 
					       double* u_grad_trial, 
					       double* u_test_dV, 
					       double* u_grad_test_dV, 
					       double* velocity,
					       double* q_m_last, 
					       double* q_cfl,
					       int* sdInfo_u_rowptr, int* sdInfo_u_colind,
					       int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
					       double* globalJacobian,
					       int nExteriorElementBoundaries_global,
					       int nQuadraturePoints_elementBoundary,
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
  using namespace LADR1Dellam;
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
	    f[nSpace],df[nSpace],a[nSpace];//todo make full tensor
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
          evaluateCoefficients_c(theta,
				 &velocity[eN_k_nSpace],
				 aL,
				 Dm,
				 u,
				 m,
				 dm,
				 f,
				 df,
				 a);

 	  for(int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_k_i=eN_k*nDOF_test_element+i;
	      int eN_k_i_nSpace=eN_k_i*nSpace;
	      for(int j=0;j<nDOF_trial_element;j++) 
		{ 
		  int eN_k_j=eN_k*nDOF_trial_element+j;
		  int eN_k_j_nSpace = eN_k_j*nSpace;
		  
		  //mwf debug
		  //std::cout<<"ellam jack eN= "<<eN<<" k= "<<k<<" dm= "<<dm<<" dt= "<<x_dt[eN_k]<<std::endl;
		  elementJacobian_u_u[i][j] += MassJacobian_weak_c(dm,u_trial[eN_k_j],u_test_dV[eN_k_i]) + 
		    SimpleDiffusionJacobian_weak_c(x_dt[eN_k],sdInfo_u_rowptr,sdInfo_u_colind,a,&u_grad_trial[eN_k_j_nSpace],
						   &u_grad_test_dV[eN_k_i_nSpace]); 		  
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
      register int eN  = elementBoundaryElementsArray[ebN*2+0];
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
	    a_ext[nSpace],//make full tensor
	    dflux_u_u_ext=0.0,
	    bc_u_ext=0.0,
	    bc_grad_u_ext[nSpace],
	    bc_m_ext=0.0,
	    bc_dm_ext=0.0,
	    bc_f_ext[nSpace],
	    bc_df_ext[nSpace],
	    bc_a_ext[nSpace],//make full tensor
	    fluxJacobian_u_u[nDOF_trial_element];
	  // 
	  //calculate the solution and gradients at quadrature points 
	  // 
	  u_ext=0.0;
	  for (int I=0;I<nSpace;I++)
	    {
	      grad_u_ext[I] = 0.0;
	      bc_grad_u_ext[I] = 0.0;
	    }
	  for (int j=0;j<nDOF_trial_element;j++) 
	    { 
	      register int eN_j = eN*nDOF_trial_element+j,
		ebNE_kb_j = ebNE_kb*nDOF_trial_element+j,
		ebNE_kb_j_nSpace= ebNE_kb_j*nSpace;
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
	  //calculate the internal and external trace of the pde coefficients 
	  // 
	  //todo add aT make full tensor
	  evaluateCoefficients_c(theta,
				 &ebqe_velocity_ext[ebNE_kb_nSpace],
				 aL,
				 Dm,
				 u_ext,
				 m_ext,
				 dm_ext,
				 f_ext,
				 df_ext,
				 a_ext);
	  evaluateCoefficients_c(theta,
				 &ebqe_velocity_ext[ebNE_kb_nSpace],
				 aL,
				 Dm,
				 bc_u_ext,
				 bc_m_ext,
				 bc_dm_ext,
				 bc_f_ext,
				 bc_df_ext,
				 bc_a_ext);
	  // 
	  //calculate the numerical fluxes 
	  // 
	  exteriorOutflowFluxDerivative_c(&ebqe_n_ext[ebNE_kb_nSpace],
					  &ebqe_velocity_ext[ebNE_kb_nSpace],
					  dflux_u_u_ext);
	  //
	  //make sure have correct dt scaling for boundaries
	  //trapezoidal rule so use 0.5dt
	  dflux_u_u_ext *= 0.5*dtnp1;

	  //calculate the flux jacobian
	  //
	  for (int j=0;j<nDOF_trial_element;j++)
	    {
	      register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
	      
	      fluxJacobian_u_u[j]=ExteriorNumericalAdvectiveFluxJacobian_c(dflux_u_u_ext,u_trial_ext[ebNE_kb_j]);
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
		  //mwf debug
		  //std::cout<<"LADR1Dellam eN= "<<eN<<" ebN= "<<ebN<<" fluxJacobian_u_u["<<j<<"]= "<<fluxJacobian_u_u[j]<<std::endl;
		}//j
	    }//i
	}//kb
    }//ebNE
}//computeJacobian

/**********************************************************************
  just tag points that need to be integrated forward for inflow boundary
  need to decide how to handle transient velocity field. Probably
  pass in new time level and old time level one next 
  Right now ignores DOF or Flux boundary flag, selects for tracking if
  velocity is inflow
 **********************************************************************/
extern "C" void markInflowBoundaryPoints(double tn, 
					 double tnp1,
					 double t,
					 int nExteriorElementBoundaries_global,
					 int nQuadraturePoints_elementBoundary,
					 const int* exteriorElementBoundariesArray,
					 const int* elementBoundaryElementsArray,
					 const int* elementBoundaryLocalElementBoundariesArray,
					 const double* ebqe_x,
					 const double* ebqe_n_ext,
					 const double* ebqe_velocity_ext_last,
					 const double* ebqe_velocity_ext,
					 const int* isDOFBoundary_u,
					 const int* isFluxBoundary_u,
					 int* element_track,//which element is point in
					 int* flag_track)  //>=0 track, -1 don't

{
  using namespace LADR1Dellam;
  for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    {
      const int ebN = exteriorElementBoundariesArray[ebNE];
      const int eN  = elementBoundaryElementsArray[ebN*2+0];
      for (int kb=0; kb < nQuadraturePoints_elementBoundary; kb++)
	{
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
	    ebNE_kb_nSpace = ebNE_kb*nSpace;
	  int needToTrackPoint=0;
	  tagInflowPointForTracking_c(tn,tnp1,t,
				      &ebqe_n_ext[ebNE_kb_nSpace],
				      &ebqe_velocity_ext_last[ebNE_kb_nSpace],
				      &ebqe_velocity_ext[ebNE_kb_nSpace],
				      needToTrackPoint);
	  flag_track[ebNE_kb] = needToTrackPoint;
	  element_track[ebNE_kb]= eN;
	}
    }

}

/**********************************************************************
  go through inflow boundary points that have been tracked forward
  and accumulate the inflow flux into the residual for the test functions
  with non-zero support in those areas
 **********************************************************************/
extern "C" void accumulateInflowFluxInGlobalResidual(int nElements_global,//mesh representation
						     int nNodes_global,
						     int nNodes_element,
						     int nElementBoundaries_element,
						     int nExteriorElementBoundaries_global,
						     int nQuadraturePoints_elementBoundary,
						     const double * nodeArray,
						     const int * elementNodesArray,
						     const int * elementNeighborsArray, //local boundary id is associated with node across from boundary 
						     const int* exteriorElementBoundariesArray,
						     const int* elementBoundaryElementsArray,
						     const int* elementBoundaryLocalElementBoundariesArray,
						     double tp, //time level for integration points
						     double timeWeight, //temporal integration weight (uniform for now)
						     const double* dS,  //spatial integration weights
						     const double* x_track_ext, //location of forward tracked integration points, assumed for now 
                                                                                //nExteriorElementBoundaries_global x nQuadraturePoints_elementBoundary
						     const double* t_track_ext,     //time forward tracked points stopped (t^{n+1} or earlier if exited domain)
						     const int* element_track_ext,  //element each forward tracked point ended up in
						     const int* flag_track_ext,     //id for each point, 0 -- interior, 1 exited domain, -1 didn't track for some reason
						     const int* u_l2g, 
						     const double* u_dof,
						     double* q_elementResidual_u, 
						     const int* sdInfo_u_rowptr, 
						     const int* sdInfo_u_colind,
						     int offset_u, int stride_u, 
						     const int* isFluxBoundary_u,
						     const double* ebqe_bc_flux_u_ext,
						     double* globalResidual)
{

  using namespace LADR1Dellam;
  /***********************************************************************
    for current boundary integration point (x_k,t_p)
       1. get boundary flux, \sigma^{p}_k
       2. get spatial integration weight W_k, temporal weight \Delta t^p
       3. get forward tracked integration point x^{n+1}_k
       4. evaluate test functions with non-zero support at x^{n+1}_k
           a. get element eN^{n+1}_k from element_track_ext[k]
           b. use physical space definition in terms of barycentric coords,
               \lambda_0 = (x-x_1)/(x_0-x_1), \lambda_1 = 1-\lambda_0
              and assume w_i = \lambda_i on eN^{n+1}_k, i=0,...nDOF_test_element-1
       5. accumulate \sigma^{p}_k*\Delta t^p* W_k*w_i(x^{n+1}_{k}) to global residual I(i) 
  **********************************************************************/
  //todo switch to a loop over boundary integration points?
  for(int ebNE=0;ebNE<nExteriorElementBoundaries_global;ebNE++)
    {
      register int ebN = exteriorElementBoundariesArray[ebNE], 
	eN  = elementBoundaryElementsArray[ebN*2+0];
      
      //declare local storage for element residual and initialize
      register double w[nDOF_test_element];
      //loop over quadrature points, may end up on elements with different support
      for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
        {
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      w[i] = 0.0;
	    }
	  //compute indeces and declare local storage
	  register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb;
	  if (flag_track_ext[ebNE_kb] >= 0 && isFluxBoundary_u[ebNE_kb])
	    {
	      //sigma_k,W_k
	      register double totalFlux = ebqe_bc_flux_u_ext[ebNE_kb],
		spatialWeight = dS[ebNE_kb];
	      register int eN_track = element_track_ext[ebNE_kb];
	      //mwf debug
	      //std::cout<<"ellam tracked boundary point ["<<ebNE<<","<<kb<<"] --> "<<x_track[ebNE_kb*3]
	      //	       <<" eN_track= "<<eN_track<<" x0= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3]
	      //       <<" x1= "<<nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3]<<std::endl;
	      evaluateTestFunctionsOnElement(&x_track_ext[ebNE_kb*3],
					     &nodeArray[elementNodesArray[eN_track*nNodes_element+0]*3],
					     &nodeArray[elementNodesArray[eN_track*nNodes_element+1]*3],
					     w);
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int eN_track_i=eN_track*nDOF_test_element+i;
		  q_elementResidual_u[eN_track_i] += totalFlux*w[i]*spatialWeight*timeWeight;
		  globalResidual[offset_u+stride_u*u_l2g[eN_track_i]] += totalFlux*w[i]*spatialWeight*timeWeight; 
		}//i
	      
	    }//needed to track point
	    
	}//integration point per element boundary
    }//element boundaries
  
}

