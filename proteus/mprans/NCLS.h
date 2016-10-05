#ifndef NCLS_H
#define NCLS_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"

//ENTROPY FUNCTIONS and SOME FLAGS (MQL)//
#define entropy_power 2. // HERE uL and uR are dummy variables
#define ENTROPY(phi) 1./entropy_power*std::pow(std::abs(phi),entropy_power)
#define ENTROPY_PRIME(phi) std::pow(std::abs(phi),entropy_power-1.)
#define ENTROPY_GRAD(phi,phix) std::pow(std::abs(phi),entropy_power-1.)*phix*(phi>=0 ? 1. : -1.)
/////////////////
// LOG ENTROPY //
/////////////////
#define ENTROPY_LOG(phi,phiL,phiR) std::log(std::abs((phi-phiL)*(phiR-phi))+1E-14)
#define ENTROPY_LOG_PRIME(phi,phiL,phiR) (phiL+phiR-2*phi)*((phi-phiL)*(phiR-phi)>0 ? 1. : -1.)*((phi-phiL)*(phiR-phi)==0 ? 0. : 1.)/(std::abs((phi-phiL)*(phiR-phi))+1E-14) 
#define ENTROPY_LOG_GRAD(phi,phix,phiL,phiR) (phiL+phiR-2*phi)*phix*((phi-phiL)*(phiR-phi)>0 ? 1 : -1)*((phi-phiL)*(phiR-phi)==0 ? 0 : 1)/(std::abs((phi-phiL)*(phiR-phi))+1E-14) 

#define tol_sign 0.1
 
namespace proteus
{
  class NCLS_base
  {
    //The base class defining the interface
  public:
    virtual ~NCLS_base(){}
    virtual void calculateResidual(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   double* meshVelocity_dof,
				   double MOVING_DOMAIN,
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
			           double useMetrics, 
				   double alphaBDF,
				   int lag_shockCapturing, /*mwf not used yet*/
				   double shockCapturingDiffusion,
		                   double sc_uref, double sc_alpha,
				   int* u_l2g, 
				   double* elementDiameter,
				   double* u_dof,double* u_dof_old,	
				   double* velocity,
				   double* q_m,
				   double* q_u,				   
				   double* q_n,
				   double* q_dH,
				   double* q_m_betaBDF,
                                   double* q_dV,
                                   double* q_dV_last,
				   double* cfl,
				   double* q_numDiff_u, 
				   double* q_numDiff_u_last, 
				   int offset_u, int stride_u, 
				   double* globalResidual,
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   double* ebqe_velocity_ext,
				   int* isDOFBoundary_u,
				   double* ebqe_rd_u_ext,				   
				   double* ebqe_bc_u_ext,				   
				   double* ebqe_u,
				   // PARAMETERS FOR ENTROPY_VISCOSITY 
				   double cE,
				   double cMax, 
				   int ENTROPY_VISCOSITY, 
				   int IMPLICIT, 
				   int SUPG,
				   // PARAMETERS FOR LS-COUPEZ 
				   int LS_COUPEZ,
				   // PARAMETERS FOR LOG BASED ENTROPY FUNCION
				   double uL, 
				   double uR)=0;
    virtual void calculateJacobian(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   double* mesh_velocity_dof,
				   double MOVING_DOMAIN,
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
			           double useMetrics, 
				   double alphaBDF,
				   int lag_shockCapturing,/*mwf not used yet*/
				   double shockCapturingDiffusion,
				   int* u_l2g,
				   double* elementDiameter,
				   double* u_dof, 
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
				   double* ebqe_velocity_ext,
				   int* isDOFBoundary_u,
				   double* ebqe_rd_u_ext,
				   double* ebqe_bc_u_ext,
				   int* csrColumnOffsets_eb_u_u,
				   int IMPLICIT, 
				   int SUPG)=0;
    virtual void calculateWaterline(//element
                                   int* wlc,
	                           double* waterline,
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   double* meshVelocity_dof,
				   double MOVING_DOMAIN,
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
			           double useMetrics, 
				   double alphaBDF,
				   int lag_shockCapturing, /*mwf not used yet*/
				   double shockCapturingDiffusion,
		                   double sc_uref, double sc_alpha,
				   int* u_l2g, 
				   double* elementDiameter,
				   double* u_dof,double* u_dof_old,	
				   double* velocity,
				   double* q_m,
				   double* q_u,				   
				   double* q_n,
				   double* q_dH,
				   double* q_m_betaBDF,
				   double* cfl,
				   double* q_numDiff_u, 
				   double* q_numDiff_u_last, 
				   int offset_u, int stride_u, 
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   int* elementBoundaryMaterialTypes,
				   double* ebqe_velocity_ext,
				   int* isDOFBoundary_u,
				   double* ebqe_bc_u_ext,				   
				   double* ebqe_u)=0;	
  };
  
  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class NCLS : public NCLS_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
  NCLS():      
    nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
    ck()
    {}
    
    inline void evaluateCoefficients(const double v[nSpace],
				     const double& u,
				     const double grad_u[nSpace],
				     double& m,
				     double& dm,
				     double& H,
				     double dH[nSpace])
    {
      m = u;
      dm=1.0;
      H = 0.0;
      for (int I=0; I < nSpace; I++)
	{
	  H += v[I]*grad_u[I];
	  dH[I] = v[I];
	}
    }

    inline
    void calculateCFL(const double& elementDiameter,
		      const double dH[nSpace],
		      double& cfl)
    {
      double h,nrm_v;
      h = elementDiameter;
      nrm_v=0.0;
      for(int I=0;I<nSpace;I++)
	nrm_v+=dH[I]*dH[I];
      nrm_v = sqrt(nrm_v);
      cfl = nrm_v/h;
    }

    inline
    void calculateSubgridError_tau(const double& elementDiameter,
				   const double& dmt,
				   const double dH[nSpace],
				   double& cfl,
				   double& tau)
    {
      double h,nrm_v,oneByAbsdt;
      h = elementDiameter;
      nrm_v=0.0;
      for(int I=0;I<nSpace;I++)
	nrm_v+=dH[I]*dH[I];
      nrm_v = sqrt(nrm_v);
      cfl = nrm_v/h;
      oneByAbsdt =  fabs(dmt);
      tau = 1.0/(2.0*nrm_v/h + oneByAbsdt + 1.0e-8);
    }

 
    inline
    void calculateSubgridError_tau(     const double&  Ct_sge,
                                        const double   G[nSpace*nSpace],
					const double&  A0,
					const double   Ai[nSpace],
					double& tau_v,
					double& cfl)	
    {
      double v_d_Gv=0.0; 
      for(int I=0;I<nSpace;I++) 
         for (int J=0;J<nSpace;J++) 
           v_d_Gv += Ai[I]*G[I*nSpace+J]*Ai[J];     
    
      tau_v = 1.0/sqrt(Ct_sge*A0*A0 + v_d_Gv + 1.0e-8);    
    } 
 
    void exteriorNumericalFlux(const double n[nSpace],
			       const double& bc_u,
			       const double& u,
			       const double velocity[nSpace],
			       const double velocity_movingDomain[nSpace],
			       double& flux)
    {
      double flow_total=0.0,flow_fluid=0.0,flow_movingDomain=0.0;
      for (int I=0; I < nSpace; I++)
	{
	  flow_fluid += n[I]*velocity[I];
	  //flow_movingDomain -= n[I]*velocity_movingDomain[I];
	}
      flow_total = flow_fluid+flow_movingDomain;
      if (flow_total > 0.0)
	{
	  flux = u*flow_movingDomain;
	}
      else
	{
	  flux = bc_u*flow_movingDomain - flow_fluid*(u-bc_u);
	  //std::cout<<"bc_u "<<bc_u<<" flow_fluid "<<flow_fluid<<" u "<<u<<std::endl;
	}
    }
    
    inline
    void exteriorNumericalFluxDerivative(const double n[nSpace],
					 const double velocity[nSpace],
					 const double velocity_movingDomain[nSpace],
					 double& dflux)
    {
      double flow_total=0.0,flow_fluid=0.0,flow_movingDomain=0.0;
      for (int I=0; I < nSpace; I++)
	{
	  flow_fluid += n[I]*velocity[I];
	  //flow_movingDomain -= n[I]*velocity_movingDomain[I];
	}
      flow_total=flow_fluid+flow_movingDomain;
      if (flow_total > 0.0)
	{
	  dflux = flow_movingDomain;
	}
      else
	{
	  dflux = -flow_fluid;
	}
    }

    void calculateResidual(//element
			   double* mesh_trial_ref,
			   double* mesh_grad_trial_ref,
			   double* mesh_dof,
			   double* mesh_velocity_dof,
			   double MOVING_DOMAIN,
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
			   double useMetrics, 
			   double alphaBDF,
			   int lag_shockCapturing, /*mwf not used yet*/
			   double shockCapturingDiffusion,
			   double sc_uref, double sc_alpha,
			   int* u_l2g, 
			   double* elementDiameter,
			   double* u_dof,double* u_dof_old,			   
			   double* velocity,
			   double* q_m,
			   double* q_u,				   
			   double* q_n,
			   double* q_dH,
			   double* q_m_betaBDF,
                           double* q_dV,
                           double* q_dV_last,
			   double* cfl,
			   double* q_numDiff_u, 
			   double* q_numDiff_u_last, 
			   int offset_u, int stride_u, 
			   double* globalResidual,
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
			   double* ebqe_velocity_ext,
			   int* isDOFBoundary_u,
			   double* ebqe_rd_u_ext,
			   double* ebqe_bc_u_ext,
			   double* ebqe_u,
			   // PARAMETERS FOR ENTROPY_VISCOSITY 
			   double cE,
			   double cMax, 
			   int ENTROPY_VISCOSITY, 
			   int IMPLICIT, 
			   int SUPG,
			   // PARAMETERS FOR LS-COUPEZ
			   int LS_COUPEZ,
			   // PARAMETERS FOR LOG BASED ENTROPY FUNCTION
			   double uL, 
			   double uR)
    {
      // ** COMPUTE QUANTITIES PER CELL (MQL) ** //
      double entropy_max=-1.E10, entropy_min=1.E10, cell_entropy_mean, cell_volume, volume=0, entropy_mean=0;
      double cell_vel_max, cell_entropy_residual;
      double dt = 1./alphaBDF; // HACKED to work just for BDF1
      double entropy_residual[nElements_global], vel_max[nElements_global];
      double entropy_normalization_factor=1.0;

      if (ENTROPY_VISCOSITY==1)
	{
	  for(int eN=0;eN<nElements_global;eN++)
	    {      
	      // FOR LS-COUPEZ
	      double llambda = 0.1; //elementDiameter[eN]/dt; // h/dt
	      cell_volume = 0;
	      cell_entropy_mean = 0;
	      cell_vel_max = 0;
	      cell_entropy_residual = 0;
	      //loop over quadrature points and compute integrands
	      for  (int k=0;k<nQuadraturePoints_element;k++)
		{
		  //compute indeces and declare local storage
		  register int eN_k = eN*nQuadraturePoints_element+k,
		    eN_k_nSpace = eN_k*nSpace,
		    eN_nDOF_trial_element = eN*nDOF_trial_element;
		  register double un=0.0,unm1=0, grad_un[nSpace], vn[nSpace],
		    jac[nSpace*nSpace],jacDet,jacInv[nSpace*nSpace],
		    u_grad_trial[nDOF_trial_element*nSpace],
		    dV,x,y,z;
		  ck.calculateMapping_element(eN,k,mesh_dof,mesh_l2g,mesh_trial_ref,mesh_grad_trial_ref,jac,jacDet,jacInv,x,y,z);
		  //get the physical integration weight
		  dV = fabs(jacDet)*dV_ref[k];
		  //get the trial function gradients
		  ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
		  //get the solution at quad point at tn and tnm1
		  ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],un);
		  ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],unm1);
		  //get the solution gradients at tn
		  ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_un);
		  //velocity at tn
		  vn[0] = velocity[eN_k_nSpace];
		  vn[1] = velocity[eN_k_nSpace+1];
		  //ADJUST VELOCITY BASED ON LS-COUPEZ
		  double grad_un_norm = 0;
		  for (int I=0; I < nSpace; I++)
		    grad_un_norm += std::pow(grad_un[I],2);
		  grad_un_norm = std::sqrt(grad_un_norm);
		  //double sign_un = (un > 0 ? 1. : -1.) * (std::abs(un)==0 ? 0. : 1.);
		  double sign_un = (un > tol_sign*uR ? 1. : -1.) * (un > -tol_sign*uR ? 0. : 1.);
		  vn[0] += LS_COUPEZ*llambda*sign_un*grad_un[0]/(grad_un_norm+1E-10);
		  vn[1] += LS_COUPEZ*llambda*sign_un*grad_un[1]/(grad_un_norm+1E-10);
		  // compute vel, entropy min, max, mean and volume
		  if (LS_COUPEZ==1)
		    {
		      entropy_max = std::max(entropy_max,ENTROPY_LOG(un,uL,uR));
		      entropy_min = std::min(entropy_min,ENTROPY_LOG(un,uL,uR));
		      cell_entropy_mean += ENTROPY_LOG(un,uL,uR)*dV;
		    }
		  else
		    {
		      entropy_max = std::max(entropy_max,ENTROPY(un));
		      entropy_min = std::min(entropy_min,ENTROPY(un));
		      cell_entropy_mean += ENTROPY(un)*dV;
		    }
		  cell_volume += dV;
		  cell_vel_max = std::max(cell_vel_max,std::max(std::abs(vn[0]),std::abs(vn[1])));
		  // compute entropy residual
		  if (LS_COUPEZ==1)
		    cell_entropy_residual 
		      = std::max(std::abs(
					  (ENTROPY_LOG(un,uL,uR) - ENTROPY_LOG(unm1,uL,uR))/dt // TIME DERIVATIVE
					  +vn[0]*ENTROPY_LOG_GRAD(un,grad_un[0],uL,uR) // TRANSPORT TERM
					  +vn[1]*ENTROPY_LOG_GRAD(un,grad_un[1],uL,uR)
					  -LS_COUPEZ*llambda*sign_un*ENTROPY_LOG_PRIME(un,uL,uR)*(1-std::pow(un/uR,2)) // FORCE TERM IN LS-COUPEZ
					  )
				 ,cell_entropy_residual);
		  else
		    cell_entropy_residual 
		      = std::max(std::abs(
					  (ENTROPY(un) - ENTROPY(unm1))/dt // TIME DERIVATIVE
					  +vn[0]*ENTROPY_GRAD(un,grad_un[0]) // TRANSPORT TERM
					  +vn[1]*ENTROPY_GRAD(un,grad_un[1])
					  -LS_COUPEZ*llambda*sign_un*ENTROPY_PRIME(un)*(1-std::pow(un/uR,2)) // FORCE TERM IN LS-COUPEZ
					  )
				 ,cell_entropy_residual);
		}
	      volume += cell_volume;
	      entropy_mean += cell_entropy_mean;
	      vel_max[eN]=cell_vel_max;
	      entropy_residual[eN] = cell_entropy_residual;
	    }//elements
	  entropy_mean /= volume;
	  // ** END OF CELL COMPUTATIONS (MQL) ** //
	  entropy_normalization_factor = std::max(std::abs(entropy_max-entropy_mean),
						  std::abs(entropy_min-entropy_mean));
	} //ENTROPY_VISCOSITY==1
      //cek should this be read in?
      double Ct_sge = 4.0;

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
	  // FOR LS-COUPEZ
	  double llambda = 0.1; //elementDiameter[eN]/dt; // h/dt
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
	      register double u=0.0, u_old=0.0, u_star=0, grad_u[nSpace], grad_u_old[nSpace], grad_u_star[nSpace],
		vn[nSpace],
		m_star=0.0, dm_star=0.0, m=0.0,dm=0.0,
		H_star=0.0, H=0.0, dH_star[nSpace], dH[nSpace],
		f_star[nSpace],df_star[nSpace],//for MOVING_DOMAIN
		m_t=0.0,dm_t=0.0,
		pdeResidual_u_star=0.0,
		Lstar_u[nDOF_test_element],
		subgridError_u_star=0.0,
		tau=0.0,tau0=0.0,tau1=0.0,
		numDiff0=0.0,numDiff1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		u_grad_trial[nDOF_trial_element*nSpace],
		u_test_dV[nDOF_trial_element],
		u_grad_test_dV[nDOF_test_element*nSpace],
		dV,x,y,z,xt,yt,zt,
		G[nSpace*nSpace],G_dd_G,tr_G;//,norm_Rv;
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
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.calculateG(jacInv,G,G_dd_G,tr_G);
	      //get the trial function gradients
	      ck.gradTrialFromRef(&u_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,u_grad_trial);
	      //get the solution
	      ck.valFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u);
	      ck.valFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],&u_trial_ref[k*nDOF_trial_element],u_old);
	      //get the solution gradients
	      ck.gradFromDOF(u_dof,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u);
	      ck.gradFromDOF(u_dof_old,&u_l2g[eN_nDOF_trial_element],u_grad_trial,grad_u_old);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  u_test_dV[j] = u_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      u_grad_test_dV[j*nSpace+I]   = u_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		    }
		}
	      //save solution at quadrature points for other models to use
	      q_u[eN_k]=u;
	      for (int I=0;I<nSpace;I++)
		q_n[eN_k_nSpace+I]  = grad_u[I];

	      // COMPUTE u and u_grad star to allow easy change between BACKWARD OR FORWARD EULER (for transport)
	      u_star = IMPLICIT*u+(1-IMPLICIT)*u_old;
	      for (int I=0; I<nSpace; I++)
		grad_u_star[I] = IMPLICIT*grad_u[I]+(1-IMPLICIT)*grad_u_old[I];

	      //
	      //calculate pde coefficients at quadrature points
	      //
	      // ADJUST VELOCITY FOR LS-COUPEZ
	      double grad_u_star_norm = 0;
	      for (int I=0; I < nSpace; I++)
		grad_u_star_norm += std::pow(grad_u_star[I],2);
	      grad_u_star_norm = std::sqrt(grad_u_star_norm);
	      //double sign_u_star = (u_star > 0 ? 1. : -1.) * (std::abs(u_star)==0 ? 0. : 1.);
	      double sign_u_star = (u_star > tol_sign*uR ? 1. : -1.) * (u_star > -tol_sign*uR ? 0. : 1.);
	      vn[0] = velocity[eN_k_nSpace];
	      vn[1] = velocity[eN_k_nSpace+1];
	      vn[0] += LS_COUPEZ*llambda*sign_u_star*grad_u_star[0]/(grad_u_star_norm+1E-10);
	      vn[1] += LS_COUPEZ*llambda*sign_u_star*grad_u_star[1]/(grad_u_star_norm+1E-10);

	      //evaluateCoefficients(&velocity[eN_k_nSpace],
	      evaluateCoefficients(vn,
				   u,
				   grad_u,
				   m,
				   dm,
				   H,
				   dH); //we don't need H and dH but we DO need m (to compute m_t)
	      //evaluateCoefficients(&velocity[eN_k_nSpace],
	      evaluateCoefficients(vn, // computed with adjusted velocity based on LS_COUPEZ
				   u_star,
				   grad_u_star,
				   m_star,
				   dm_star,
				   H_star, // H_star = velocity*grad(u_star)
				   dH_star); // dH_star = velocity
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      for (int I=0;I<nSpace;I++)
		{
		  f_star[I] = -MOVING_DOMAIN*m_star*mesh_velocity[I];
		  df_star[I] = -MOVING_DOMAIN*dm_star*mesh_velocity[I];
		}
	      //
	      //calculate time derivative at quadrature points
	      //
              if (q_dV_last[eN_k] <= -100)
                q_dV_last[eN_k] = dV;
              q_dV[eN_k] = dV;
	      ck.bdf(alphaBDF,
		     q_m_betaBDF[eN_k],
		     m,
		     dm,
		     m_t,
		     dm_t);
	      if (ENTROPY_VISCOSITY==0)
		{
		  //
		  //calculate subgrid error (strong residual and adjoint)
		  //
		  //calculate strong residual
		  pdeResidual_u_star = ck.Mass_strong(m_t) +
		    ck.Hamiltonian_strong(dH_star,grad_u)+
		    MOVING_DOMAIN*ck.Advection_strong(df_star,grad_u);//cek I don't think all mesh motion will be divergence free so we may need to go add the divergence
	      
		  //calculate adjoint
		  for (int i=0;i<nDOF_test_element;i++)
		    {
		      //register int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
		      //Lstar_u[i]  = ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[eN_k_i_nSpace]);
		      register int i_nSpace = i*nSpace;
		      Lstar_u[i]  = ck.Hamiltonian_adjoint(dH_star,&u_grad_test_dV[i_nSpace]) 
			+ MOVING_DOMAIN*ck.Advection_adjoint(df_star,&u_grad_test_dV[i_nSpace]);
		    }
		  //calculate tau and tau*Res
		  double subgridErrorVelocity[nSpace];
		  for (int I=0;I<nSpace;I++)
		    subgridErrorVelocity[I] = dH_star[I] - MOVING_DOMAIN*df_star[I];
		  
		  calculateSubgridError_tau(elementDiameter[eN],
					    dm_t,
					    subgridErrorVelocity,//dH,
					    cfl[eN_k],
					    tau0);
		  calculateSubgridError_tau(Ct_sge,
					    G,
					    dm_t,
					    subgridErrorVelocity,//dH,
					    tau1,
					    cfl[eN_k]);
		  tau = useMetrics*tau1+(1.0-useMetrics)*tau0;
		  subgridError_u_star = -tau*pdeResidual_u_star;
		  //
		  //calculate shock capturing diffusion
		  //	      
		  ck.calculateNumericalDiffusion(shockCapturingDiffusion,elementDiameter[eN],pdeResidual_u_star,grad_u_star,numDiff0);	      
		  //ck.calculateNumericalDiffusion(shockCapturingDiffusion,G,pdeResidual_u,grad_u_old,numDiff1);
		  ck.calculateNumericalDiffusion(shockCapturingDiffusion,sc_uref, sc_alpha,G,G_dd_G,pdeResidual_u_star,grad_u_star,numDiff1);
		  q_numDiff_u[eN_k] = useMetrics*numDiff1+(1.0-useMetrics)*numDiff0*0;
		} // ENTROPY_VISCOSITY=0
	      else 
		{
		  // CALCULATE CFL //
		  calculateCFL(elementDiameter[eN],dH_star,cfl[eN_k]); // TODO: ADJUST SPEED IF MESH IS MOVING
		  // ** LINEAR DIFFUSION (MQL) ** //
		  // calculate linear viscosity 
		  double h=elementDiameter[eN];
		  double vMax = std::max(std::abs(dH_star[0]),std::abs(dH_star[1]));
		  double linear_viscosity = cMax*h*vel_max[eN]; // Cell based
		  
		  // ** ENTROPY VISCOSITY (MQL) ** //
		  // calculate entropy residual
		  double entropy_viscosity = cE*h*h*entropy_residual[eN]/entropy_normalization_factor;
		  q_numDiff_u[eN_k] = std::min(linear_viscosity,entropy_viscosity);
		  //q_numDiff_u[eN_k] = linear_viscosity;
		}
	      //update element residual 
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  //register int eN_k_i=eN_k*nDOF_test_element+i,
		   // eN_k_i_nSpace = eN_k_i*nSpace;
		  register int  i_nSpace=i*nSpace;
		  elementResidual_u[i] += 
		    dt*ck.Mass_weak(m_t,u_test_dV[i]) + 
		    dt*ck.Hamiltonian_weak(H_star,u_test_dV[i]) + 
		    dt*MOVING_DOMAIN*ck.Advection_weak(f_star,&u_grad_test_dV[i_nSpace])+
		    dt*SUPG*ck.SubgridError(subgridError_u_star,Lstar_u[i]) + 
		    dt*ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u_star,&u_grad_test_dV[i_nSpace])
		    //-dt*LS_COUPEZ*llambda*sign_u_star*(1-std::sqrt(std::min(std::pow(u_star/uR,2),1.0)))*u_test_dV[i] // FORCE TERM IN LS-COUPEZ
		    -dt*LS_COUPEZ*llambda*sign_u_star*(1-std::pow(u_star/uR,2))*u_test_dV[i] // FORCE TERM IN LS-COUPEZ
		    ;
		}//i
	      //
	      //cek/ido todo, get rid of m, since u=m
	      //save momentum for time history and velocity for subgrid error
	      //save solution for other models 
	      //
	      q_m[eN_k] = m;
	      q_u[eN_k] = u;
	      for (int I=0;I<nSpace;I++)
		{
		  int eN_k_I = eN_k*nSpace+I;
		  q_dH[eN_k_I] = dH[I];                     
		}
	    }
	  //
	  //load element into global residual and save element residual
	  //
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_i=eN*nDOF_test_element+i;
          
	      globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];
	    }//i
	}//elements
      //
      //loop over exterior element boundaries to calculate surface integrals and load into element and global residuals
      //
      //ebNE is the Exterior element boundary INdex
      //ebN is the element boundary INdex
      //eN is the element index
      //std::cout <<nExteriorElementBoundaries_global<<std::endl;
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
	  register int ebN = exteriorElementBoundariesArray[ebNE], 
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_u[nDOF_test_element];
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
		H_ext=0.0,
		dH_ext[nSpace],
		//f_ext[nSpace],//MOVING_DOMAIN
		//df_ext[nSpace],//MOVING_DOMAIN
		//flux_ext=0.0,
		bc_u_ext=0.0,
		bc_grad_u_ext[nSpace],
		bc_m_ext=0.0,
		bc_dm_ext=0.0,
		bc_H_ext=0.0,
		bc_dH_ext[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,
		u_test_dS[nDOF_test_element],
		u_grad_trial_trace[nDOF_trial_element*nSpace],
		normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		G[nSpace*nSpace],G_dd_G,tr_G,flux_ext;
	      // 
	      //calculate the solution and gradients at quadrature points 
	      // 
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
	      ck.calculateMappingVelocity_elementBoundary(eN,
							  ebN_local,
							  kb,
							  ebN_local_kb,
							  mesh_velocity_dof,
							  mesh_l2g,
							  mesh_trial_trace_ref,
							  xt_ext,yt_ext,zt_ext,
							  normal,
							  boundaryJac,
							  metricTensor,
							  integralScaling);
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
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
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*ebqe_rd_u_ext[ebNE_kb];
	      // 
	      //calculate the pde coefficients using the solution and the boundary values for the solution 
	      // 
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   u_ext,
				   grad_u_ext,
				   m_ext,
				   dm_ext,
				   H_ext,
				   dH_ext);
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   bc_u_ext,
				   bc_grad_u_ext,
				   bc_m_ext,
				   bc_dm_ext,
				   bc_H_ext,
				   bc_dH_ext);
	      //
	      //moving mesh
	      //
	      double velocity_ext[nSpace];
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      for (int I=0;I<nSpace;I++)
		velocity_ext[I] = - MOVING_DOMAIN*mesh_velocity[I];
	      // 
	      //calculate the numerical fluxes 
	      // 
	      exteriorNumericalFlux(normal,
				    bc_u_ext,
				    u_ext,
				    dH_ext,
				    velocity_ext,
				    flux_ext);
	      ebqe_u[ebNE_kb] = u_ext;
	      
	      //std::cout<<u_ext<<ebqe_bc_u_ext
	      
	      //
	      //update residuals
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  //int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;

		  elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_ext,u_test_dS[i]);  		  
		}//i
	    }//kb
	  //
	  //update the element and global residual storage
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;

	      //globalResidual[offset_u+stride_u*u_l2g[eN_i]] += MOVING_DOMAIN*elementResidual_u[i];
              globalResidual[offset_u+stride_u*u_l2g[eN_i]] += elementResidual_u[i];	  
	      
	    }//i
	}//ebNE
    }
    
    void calculateJacobian(//element
			   double* mesh_trial_ref,
			   double* mesh_grad_trial_ref,
			   double* mesh_dof,
			   double* mesh_velocity_dof,
			   double MOVING_DOMAIN,
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
			   double useMetrics, 
			   double alphaBDF,
			   int lag_shockCapturing,/*mwf not used yet*/
			   double shockCapturingDiffusion,
			   int* u_l2g,
			   double* elementDiameter,
			   double* u_dof, 
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
			   double* ebqe_velocity_ext,
			   int* isDOFBoundary_u,
			   double* ebqe_rd_u_ext,
			   double* ebqe_bc_u_ext,
			   int* csrColumnOffsets_eb_u_u,
			   int IMPLICIT, 
			   int SUPG)
    {
      double dt = 1./alphaBDF; // valid just for forward/backward euler
      double Ct_sge = 4.0;
    
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
		H=0.0,dH[nSpace],
		f[nSpace],df[nSpace],//MOVING_MESH
		m_t=0.0,dm_t=0.0,
		dpdeResidual_u_u[nDOF_trial_element],
		Lstar_u[nDOF_test_element],
		dsubgridError_u_u[nDOF_trial_element],
		tau=0.0,tau0=0.0,tau1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		u_grad_trial[nDOF_trial_element*nSpace],
		dV,
		u_test_dV[nDOF_test_element],
		u_grad_test_dV[nDOF_test_element*nSpace],
		x,y,z,xt,yt,zt,
		G[nSpace*nSpace],G_dd_G,tr_G;
	      //
	      //calculate solution and gradients at quadrature points
	      //
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
	      ck.calculateMappingVelocity_element(eN,
						  k,
						  mesh_velocity_dof,
						  mesh_l2g,
						  mesh_trial_ref,
						  xt,yt,zt);
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
	      evaluateCoefficients(&velocity[eN_k_nSpace],
				   u,
				   grad_u,
				   m,
				   dm,
				   H,
				   dH);
	      //
	      //moving mesh
	      //
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt;
	      mesh_velocity[1] = yt;
	      mesh_velocity[2] = zt;
	      for (int I=0;I<nSpace;I++)
		{
		  f[I] = -MOVING_DOMAIN*m*mesh_velocity[I];
		  df[I] = -MOVING_DOMAIN*dm*mesh_velocity[I];
		}
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
		  //int eN_k_i_nSpace = (eN_k*nDOF_trial_element+i)*nSpace;
		  //Lstar_u[i]=ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[eN_k_i_nSpace]);
		  register int i_nSpace = i*nSpace;
		  Lstar_u[i]=ck.Hamiltonian_adjoint(dH,&u_grad_test_dV[i_nSpace]) + MOVING_DOMAIN*ck.Advection_adjoint(df,&u_grad_test_dV[i_nSpace]);
	      
		}
	      //calculate the Jacobian of strong residual
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //int eN_k_j=eN_k*nDOF_trial_element+j;
		  //int eN_k_j_nSpace = eN_k_j*nSpace;
		  int j_nSpace = j*nSpace;
		  dpdeResidual_u_u[j]=ck.MassJacobian_strong(dm_t,u_trial_ref[k*nDOF_trial_element+j]) +
		    ck.HamiltonianJacobian_strong(dH,&u_grad_trial[j_nSpace]) +
		    MOVING_DOMAIN*ck.AdvectionJacobian_strong(df,&u_grad_trial[j_nSpace]);

		}
	      //tau and tau*Res
	      double subgridErrorVelocity[nSpace];
	      for (int I=0;I<nSpace;I++)
		subgridErrorVelocity[I] = dH[I] - MOVING_DOMAIN*df[I];

	      calculateSubgridError_tau(elementDiameter[eN],
					dm_t,
					subgridErrorVelocity,//dH,
	   	         		cfl[eN_k],
					tau0);
  
              calculateSubgridError_tau(Ct_sge,
                                        G,
					dm_t,
					subgridErrorVelocity,//dH,
					tau1,
				        cfl[eN_k]);
              tau = useMetrics*tau1+(1.0-useMetrics)*tau0;

	      for(int j=0;j<nDOF_trial_element;j++)
		dsubgridError_u_u[j] = -tau*dpdeResidual_u_u[j];
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  //int eN_k_i=eN_k*nDOF_test_element+i;
		  //int eN_k_i_nSpace=eN_k_i*nSpace;
		  for(int j=0;j<nDOF_trial_element;j++) 
		    { 
		      //int eN_k_j=eN_k*nDOF_trial_element+j;
		      //int eN_k_j_nSpace = eN_k_j*nSpace;
		      int j_nSpace = j*nSpace;
		      int i_nSpace = i*nSpace;
		      elementJacobian_u_u[i][j] += 
			dt*ck.MassJacobian_weak(dm_t,u_trial_ref[k*nDOF_trial_element+j],u_test_dV[i]) + 
			dt*IMPLICIT*ck.HamiltonianJacobian_weak(dH,&u_grad_trial[j_nSpace],u_test_dV[i]) +
			dt*IMPLICIT*MOVING_DOMAIN*ck.AdvectionJacobian_weak(df,u_trial_ref[k*nDOF_trial_element+j],&u_grad_test_dV[i_nSpace]) +
			dt*SUPG*IMPLICIT*ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u[i]) + 
			dt*IMPLICIT*ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&u_grad_trial[j_nSpace],&u_grad_test_dV[i_nSpace]); 
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
		H_ext=0.0,
		dH_ext[nSpace],
		bc_grad_u_ext[nSpace],
		bc_H_ext=0.0,
		bc_dH_ext[nSpace],
		//f_ext[nSpace],
		//df_ext[nSpace],
		dflux_u_u_ext=0.0,
		bc_u_ext=0.0,
		//bc_grad_u_ext[nSpace],
		bc_m_ext=0.0,
		bc_dm_ext=0.0,
		//bc_f_ext[nSpace],
		//bc_df_ext[nSpace],
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
		normal[nSpace],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
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
	      ck.calculateMappingVelocity_elementBoundary(eN,
							  ebN_local,
							  kb,
							  ebN_local_kb,
							  mesh_velocity_dof,
							  mesh_l2g,
							  mesh_trial_trace_ref,
							  xt_ext,yt_ext,zt_ext,
							  normal,
							  boundaryJac,
							  metricTensor,
							  integralScaling);
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //dS = metricTensorDetSqrt*dS_ref[kb];
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
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*ebqe_rd_u_ext[ebNE_kb];
	      // 
	      //calculate the internal and external trace of the pde coefficients 
	      // 
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   u_ext,
				   grad_u_ext,
				   m_ext,
				   dm_ext,
				   H_ext,
				   dH_ext);
	      evaluateCoefficients(&ebqe_velocity_ext[ebNE_kb_nSpace],
				   bc_u_ext,
				   bc_grad_u_ext,
				   bc_m_ext,
				   bc_dm_ext,
				   bc_H_ext,
				   bc_dH_ext);
	      //
	      //moving domain
	      //
	      double velocity_ext[nSpace];
	      double mesh_velocity[3];
	      mesh_velocity[0] = xt_ext;
	      mesh_velocity[1] = yt_ext;
	      mesh_velocity[2] = zt_ext;
	      for (int I=0;I<nSpace;I++)
		{
		  velocity_ext[I] = - MOVING_DOMAIN*mesh_velocity[I];
		}
	      // 
	      //calculate the numerical fluxes 
	      // 
	      exteriorNumericalFluxDerivative(normal,
					      dH_ext,
					      velocity_ext,
					      dflux_u_u_ext);
	      //
	      //calculate the flux jacobian
	      //
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  //register int ebNE_kb_j = ebNE_kb*nDOF_trial_element+j;
		  register int ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
	      
		  fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_u_u_ext,u_trial_trace_ref[ebN_local_kb_j]);
		}//j
	      //
	      //update the global Jacobian from the flux Jacobian
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i = eN*nDOF_test_element+i;
		  //register int ebNE_kb_i = ebNE_kb*nDOF_test_element+i;
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
                      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*u_test_dS[i];
		      //globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += MOVING_DOMAIN*fluxJacobian_u_u[j]*u_test_dS[i];
		    }//j
		}//i
	    }//kb
	}//ebNE
    }//computeJacobian

    void calculateWaterline(//element
                           int* wlc,
	                   double* waterline,
			   double* mesh_trial_ref,
			   double* mesh_grad_trial_ref,
			   double* mesh_dof,
			   double* mesh_velocity_dof,
			   double MOVING_DOMAIN,
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
			   double useMetrics, 
			   double alphaBDF,
			   int lag_shockCapturing, /*mwf not used yet*/
			   double shockCapturingDiffusion,
			   double sc_uref, double sc_alpha,
			   int* u_l2g, 
			   double* elementDiameter,
			   double* u_dof,double* u_dof_old,			   
			   double* velocity,
			   double* q_m,
			   double* q_u,				   
			   double* q_n,
			   double* q_dH,
			   double* q_m_betaBDF,
			   double* cfl,
			   double* q_numDiff_u, 
			   double* q_numDiff_u_last, 
			   int offset_u, int stride_u, 
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
                           int* elementBoundaryMaterialTypes,
			   double* ebqe_velocity_ext,
			   int* isDOFBoundary_u,
			   double* ebqe_bc_u_ext,
			   double* ebqe_u)
    {

      //  Tetrehedral elements specific extraction routine for waterline extraction
      //  Loops over boundaries and checks if boundary is infact a hull (hardwired check if mattype > 6)
      //  Extracts the nodal values of boundary triangle (4th point is dropped = hardwired assumption we are dealing with linear tet)
      //  Then computes an average value and position for both negative and positive values
      //  If both positive and negative values re found, and we are actually in a triangle containing the interface
      //  a linear interpolation of negative and positive average is reported as interface location (exact in case of linear tets) 

      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++) 
	{ 
          register int ebN = exteriorElementBoundariesArray[ebNE]; 
	  register int eN  = elementBoundaryElementsArray[ebN*2+0];
	  register int bN  = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
	    
	  if (elementBoundaryMaterialTypes[ebN] >6)
	    {
	    	    
	    double val,x,y,z;	    
	    int pos=0, neg=0;
	    double xn=0.0, yn=0.0, zn=0.0, vn=0.0; 
	    double xp=0.0, yp=0.0, zp=0.0, vp=0.0;  
	    
            for (int j=0;j<nDOF_trial_element;j++) 
	    {
	      if (j != bN) {  
	       int eN_nDOF_trial_element_j      = eN*nDOF_trial_element + j;
	       int eN_nDOF_mesh_trial_element_j = eN*nDOF_mesh_trial_element + j;
               val = u_dof[u_l2g[eN_nDOF_trial_element_j]];
	       x = mesh_dof[mesh_l2g[eN_nDOF_mesh_trial_element_j]*3+0];
	       y = mesh_dof[mesh_l2g[eN_nDOF_mesh_trial_element_j]*3+1];
	       z = mesh_dof[mesh_l2g[eN_nDOF_mesh_trial_element_j]*3+2];
	       
               if (val < 0.0)
	       {
	         neg++;
		 vn+=val;
		 xn+=x;
		 yn+=y;
		 zn+=z;
	       }
	       else
	       {
	       	 pos++;
		 vp+=val;
		 xp+=x;
		 yp+=y;
		 zp+=z;
	       }
	      }	            
            } // trail for
	    

            if ((pos > 0) && (neg > 0) )
	    {
	      vp /= pos;
	      vn /= neg;
	  
	      double alpha = vp/(vp -vn);

	      waterline[wlc[0]*3 + 0] =  alpha*(xn/neg) + (1.0-alpha)*(xp/pos);
	      waterline[wlc[0]*3 + 1] =  alpha*(yn/neg) + (1.0-alpha)*(yp/pos);
	      waterline[wlc[0]*3 + 2] =  alpha*(zn/neg) + (1.0-alpha)*(zp/pos);
	      wlc[0]++;
	      	      
	    } // end value if	 
	  
	  } // end bnd mat check
	     
	}//ebNE

        //std::cout<<"CPP WLC "<<wlc[0]<<std::endl;
    } // calcWaterline

  };//NCLS

  inline NCLS_base* newNCLS(int nSpaceIn,
			    int nQuadraturePoints_elementIn,
			    int nDOF_mesh_trial_elementIn,
			    int nDOF_trial_elementIn,
			    int nDOF_test_elementIn,
			    int nQuadraturePoints_elementBoundaryIn,
			    int CompKernelFlag)
  {
    if (nSpaceIn == 2)
      return proteus::chooseAndAllocateDiscretization2D<NCLS_base,NCLS,CompKernel>(nSpaceIn,
										   nQuadraturePoints_elementIn,
										   nDOF_mesh_trial_elementIn,
										   nDOF_trial_elementIn,
										   nDOF_test_elementIn,
										   nQuadraturePoints_elementBoundaryIn,
										   CompKernelFlag);
    else
      return proteus::chooseAndAllocateDiscretization<NCLS_base,NCLS,CompKernel>(nSpaceIn,
										 nQuadraturePoints_elementIn,
										 nDOF_mesh_trial_elementIn,
										 nDOF_trial_elementIn,
										 nDOF_test_elementIn,
										 nQuadraturePoints_elementBoundaryIn,
										 CompKernelFlag);
    /* return proteus::chooseAndAllocateDiscretization<NCLS_base,NCLS>(nSpaceIn, */
    /* 									       nQuadraturePoints_elementIn, */
    /* 									       nDOF_mesh_trial_elementIn, */
    /* 									       nDOF_trial_elementIn, */
    /* 									       nDOF_test_elementIn, */
    /* 									       nQuadraturePoints_elementBoundaryIn, */
    /* 									       CompKernelFlag); */
  }
}//proteus
#endif
