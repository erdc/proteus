#ifndef RANS3P_H
#define RANS3P_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"
#include "ModelFactory.h"
//#define COMPRESSIBLE_FORM
namespace proteus
{
  class RANS3P_base
  {
  public:
    virtual ~RANS3P_base(){}
    virtual void calculateResidual(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   double* mesh_velocity_dof,
				   double MOVING_DOMAIN,//0 or 1
				   int* mesh_l2g,
				   double* dV_ref,
				   double* p_trial_ref,
				   double* p_grad_trial_ref,
				   double* p_test_ref,
				   double* p_grad_test_ref,
				   double* vel_trial_ref,
				   double* vel_grad_trial_ref,
				   double* vel_test_ref,
				   double* vel_grad_test_ref,
				   //element boundary
				   double* mesh_trial_trace_ref,
				   double* mesh_grad_trial_trace_ref,
				   double* dS_ref,
				   double* p_trial_trace_ref,
				   double* p_grad_trial_trace_ref,
				   double* p_test_trace_ref,
				   double* p_grad_test_trace_ref,
				   double* vel_trial_trace_ref,
				   double* vel_grad_trial_trace_ref,
				   double* vel_test_trace_ref,
				   double* vel_grad_test_trace_ref,					 
				   double* normal_ref,
				   double* boundaryJac_ref,
				   //physics
				   double* elementDiameter,
				   double hFactor,
				   int nElements_global,
				   double useRBLES,
			           double useMetrics, 
				   double alphaBDF,
				   double epsFact_rho,
				   double epsFact_mu, 
				   double sigma,
				   double rho_0,
				   double nu_0,
				   double rho_1,
				   double nu_1,
				   double Ct_sge,
				   double Cd_sge,
				   double C_dc,
				   int* p_l2g, 
				   int* vel_l2g, 
				   double* C_dof,
				   double* p_dof, 
				   double* u_dof, 
				   double* v_dof, 
				   double* w_dof,
				   double* g,
				   double* q_sedimentForce,
				   //double* ebqe_sedimentForce,
				   double* phi,
				   double* normal_phi,
				   double* kappa_phi,
				   double* q_mass_acc,
				   double* q_mom_u_acc,
				   double* q_mom_v_acc,
				   double* q_mom_w_acc,
				   double* q_mass_adv,
				   double* q_mass_acc_beta_bdf,double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				   double* q_velocity_sge,
				   double* q_cfl,
				   double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
				   double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
				   int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
				   int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
				   int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
				   int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
				   int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
				   int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
				   int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
				   int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
				   int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
				   int offset_p, int offset_u, int offset_v, int offset_w, 
				   int stride_p, int stride_u, int stride_v, int stride_w, 
				   double* globalResidual,
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   double* ebqe_phi_ext,
				   double* ebqe_normal_phi_ext,
				   double* ebqe_kappa_phi_ext,
				   int* isDOFBoundary_p,
				   int* isDOFBoundary_u,
				   int* isDOFBoundary_v,
				   int* isDOFBoundary_w,
				   int* isAdvectiveFluxBoundary_p,
				   int* isAdvectiveFluxBoundary_u,
				   int* isAdvectiveFluxBoundary_v,
				   int* isAdvectiveFluxBoundary_w,
				   int* isDiffusiveFluxBoundary_u,
				   int* isDiffusiveFluxBoundary_v,
				   int* isDiffusiveFluxBoundary_w,
				   double* ebqe_bc_p_ext,
				   double* ebqe_bc_flux_mass_ext,
				   double* ebqe_bc_flux_mom_u_adv_ext,
				   double* ebqe_bc_flux_mom_v_adv_ext,
				   double* ebqe_bc_flux_mom_w_adv_ext,
				   double* ebqe_bc_u_ext,
				   double* ebqe_bc_flux_u_diff_ext,
				   double* ebqe_penalty_ext,
				   double* ebqe_bc_v_ext,
				   double* ebqe_bc_flux_v_diff_ext,
				   double* ebqe_bc_w_ext,
				   double* ebqe_bc_flux_w_diff_ext,
				   double* q_velocity,
				   double* ebqe_velocity,
				   double* flux,
				   double* elementResidual_p)=0;
    virtual void calculateJacobian(//element
				   double* mesh_trial_ref,
				   double* mesh_grad_trial_ref,
				   double* mesh_dof,
				   double* mesh_velocity_dof,
				   double MOVING_DOMAIN,
				   int* mesh_l2g,
				   double* dV_ref,
				   double* p_trial_ref,
				   double* p_grad_trial_ref,
				   double* p_test_ref,
				   double* p_grad_test_ref,
				   double* vel_trial_ref,
				   double* vel_grad_trial_ref,
				   double* vel_test_ref,
				   double* vel_grad_test_ref,
				   //element boundary
				   double* mesh_trial_trace_ref,
				   double* mesh_grad_trial_trace_ref,
				   double* dS_ref,
				   double* p_trial_trace_ref,
				   double* p_grad_trial_trace_ref,
				   double* p_test_trace_ref,
				   double* p_grad_test_trace_ref,
				   double* vel_trial_trace_ref,
				   double* vel_grad_trial_trace_ref,
				   double* vel_test_trace_ref,
				   double* vel_grad_test_trace_ref,					 
				   double* normal_ref,
				   double* boundaryJac_ref,
				   //physics
				   double* elementDiameter,
				   double hFactor,
				   int nElements_global,
				   double useRBLES,
			           double useMetrics, 
				   double alphaBDF,
				   double epsFact_rho,
				   double epsFact_mu,
				   double sigma,
				   double rho_0,
				   double nu_0,
				   double rho_1,
				   double nu_1,
				   double Ct_sge,
				   double Cd_sge,
				   double C_dg,
				   int* p_l2g, 
				   int* vel_l2g,
				   double* C_dof,double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				   double* g,
				   double* q_sedimentForce,
				   //double* ebqe_sedimentForce,
				   double* phi,
				   double* normal_phi,
				   double* kappa_phi,
				   double* q_mass_acc_beta_bdf,double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				   double* q_velocity_sge,
				   double* q_cfl,
				   double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
				   int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
				   int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
				   int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
				   int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
				   int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
				   int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
				   int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
				   int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
				   int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
				   int* csrRowIndeces_p_p,int* csrColumnOffsets_p_p,
				   int* csrRowIndeces_p_u,int* csrColumnOffsets_p_u,
				   int* csrRowIndeces_p_v,int* csrColumnOffsets_p_v,
				   int* csrRowIndeces_p_w,int* csrColumnOffsets_p_w,
				   int* csrRowIndeces_u_p,int* csrColumnOffsets_u_p,
				   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
				   int* csrRowIndeces_u_v,int* csrColumnOffsets_u_v,
				   int* csrRowIndeces_u_w,int* csrColumnOffsets_u_w,
				   int* csrRowIndeces_v_p,int* csrColumnOffsets_v_p,
				   int* csrRowIndeces_v_u,int* csrColumnOffsets_v_u,
				   int* csrRowIndeces_v_v,int* csrColumnOffsets_v_v,
				   int* csrRowIndeces_v_w,int* csrColumnOffsets_v_w,
				   int* csrRowIndeces_w_p,int* csrColumnOffsets_w_p,
				   int* csrRowIndeces_w_u,int* csrColumnOffsets_w_u,
				   int* csrRowIndeces_w_v,int* csrColumnOffsets_w_v,
				   int* csrRowIndeces_w_w,int* csrColumnOffsets_w_w,
				   double* globalJacobian,
				   int nExteriorElementBoundaries_global,
				   int* exteriorElementBoundariesArray,
				   int* elementBoundaryElementsArray,
				   int* elementBoundaryLocalElementBoundariesArray,
				   double* ebqe_phi_ext,
				   double* ebqe_normal_phi_ext,
				   double* ebqe_kappa_phi_ext,
				   int* isDOFBoundary_p,
				   int* isDOFBoundary_u,
				   int* isDOFBoundary_v,
				   int* isDOFBoundary_w,
				   int* isAdvectiveFluxBoundary_p,
				   int* isAdvectiveFluxBoundary_u,
				   int* isAdvectiveFluxBoundary_v,
				   int* isAdvectiveFluxBoundary_w,
				   int* isDiffusiveFluxBoundary_u,
				   int* isDiffusiveFluxBoundary_v,
				   int* isDiffusiveFluxBoundary_w,
				   double* ebqe_bc_p_ext,
				   double* ebqe_bc_flux_mass_ext,
				   double* ebqe_bc_flux_mom_u_adv_ext,
				   double* ebqe_bc_flux_mom_v_adv_ext,
				   double* ebqe_bc_flux_mom_w_adv_ext,
				   double* ebqe_bc_u_ext,
				   double* ebqe_bc_flux_u_diff_ext,
				   double* ebqe_penalty_ext,
				   double* ebqe_bc_v_ext,
				   double* ebqe_bc_flux_v_diff_ext,
				   double* ebqe_bc_w_ext,
				   double* ebqe_bc_flux_w_diff_ext,
				   int* csrColumnOffsets_eb_p_p,
				   int* csrColumnOffsets_eb_p_u,
				   int* csrColumnOffsets_eb_p_v,
				   int* csrColumnOffsets_eb_p_w,
				   int* csrColumnOffsets_eb_u_p,
				   int* csrColumnOffsets_eb_u_u,
				   int* csrColumnOffsets_eb_u_v,
				   int* csrColumnOffsets_eb_u_w,
				   int* csrColumnOffsets_eb_v_p,
				   int* csrColumnOffsets_eb_v_u,
				   int* csrColumnOffsets_eb_v_v,
				   int* csrColumnOffsets_eb_v_w,
				   int* csrColumnOffsets_eb_w_p,
				   int* csrColumnOffsets_eb_w_u,
				   int* csrColumnOffsets_eb_w_v,
				   int* csrColumnOffsets_eb_w_w)=0;
    virtual void calculateVelocityAverage(int nExteriorElementBoundaries_global,
    					  int* exteriorElementBoundariesArray,
    					  int nInteriorElementBoundaries_global,
    					  int* interiorElementBoundariesArray,
    					  int* elementBoundaryElementsArray,
    					  int* elementBoundaryLocalElementBoundariesArray,
    					  double* mesh_dof,
    					  int* mesh_l2g,
    					  double* mesh_trial_trace_ref,
    					  double* mesh_grad_trial_trace_ref,
    					  double* normal_ref,
    					  double* boundaryJac_ref,
    					  int* vel_l2g,
    					  double* u_dof,
    					  double* v_dof,
    					  double* w_dof,
    					  double* vel_trial_trace_ref,
    					  double* ebqe_velocity,
    					  double* velocityAverage)=0;
  };
  
  template<class CompKernelType,
	   int nSpace,
	   int nQuadraturePoints_element,
	   int nDOF_mesh_trial_element,
	   int nDOF_trial_element,
	   int nDOF_test_element,
	   int nQuadraturePoints_elementBoundary>
  class RANS3P : public RANS3P_base
  {
  public:
    const int nDOF_test_X_trial_element;
    CompKernelType ck;
    RANS3P():
      nDOF_test_X_trial_element(nDOF_test_element*nDOF_trial_element),
      ck()
    {}

    inline double smoothedHeaviside(double eps, double phi)
    {
      double H;
      if (phi > eps)
	H=1.0;
      else if (phi < -eps)
	H=0.0;
      else if (phi==0.0)
	H=0.5;
      else
	H = 0.5*(1.0 + phi/eps + sin(M_PI*phi/eps)/M_PI);
      return H;
    }

    inline double smoothedHeaviside_integral(double eps, double phi)
    {
      double HI;
      if (phi > eps)
	{
	  HI= phi - eps							\
	    + 0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) \
	    - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
	}
      else if (phi < -eps)
	{
	  HI=0.0;
	}
      else
	{
	  HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) \
	    - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
	}
      return HI;
    }
 
    inline double smoothedDirac(double eps, double phi)
    {
      double d;
      if (phi > eps)
	d=0.0;
      else if (phi < -eps)
	d=0.0;
      else
	d = 0.5*(1.0 + cos(M_PI*phi/eps))/eps;
      return d;
    }

    inline
    void evaluateCoefficients(const double eps_rho,
			      const double eps_mu,
			      const double sigma,
			      const double rho_0,
			      const double nu_0,
			      const double rho_1,
			      const double nu_1,
			      const double g[nSpace],
			      const double sedimentForce[nSpace],
			      const double& phi,
			      const double n[nSpace],
			      const double& kappa,
			      const double& C,
			      const double grad_C[nSpace],
			      const double& p,
			      const double grad_p[nSpace],
			      const double& u,
			      const double& v,
			      const double& w,
			      double& mass_acc,
			      double& mom_u_acc,
			      double& dmom_u_acc_u,
			      double& mom_v_acc,
			      double& dmom_v_acc_v,
			      double& mom_w_acc,
			      double& dmom_w_acc_w,
			      double mass_adv[nSpace],
			      double dmass_adv_u[nSpace],
			      double dmass_adv_v[nSpace],
			      double dmass_adv_w[nSpace],
			      double mom_u_adv[nSpace],
			      double dmom_u_adv_u[nSpace],
			      double dmom_u_adv_v[nSpace],
			      double dmom_u_adv_w[nSpace],
			      double mom_v_adv[nSpace],
			      double dmom_v_adv_u[nSpace],
			      double dmom_v_adv_v[nSpace],
			      double dmom_v_adv_w[nSpace],
			      double mom_w_adv[nSpace],
			      double dmom_w_adv_u[nSpace],
			      double dmom_w_adv_v[nSpace],
			      double dmom_w_adv_w[nSpace],
			      double mom_u_diff_ten[nSpace],
			      double mom_v_diff_ten[nSpace],
			      double mom_w_diff_ten[nSpace],
			      double mom_uv_diff_ten[1],
			      double mom_uw_diff_ten[1],
			      double mom_vu_diff_ten[1],
			      double mom_vw_diff_ten[1],
			      double mom_wu_diff_ten[1],
			      double mom_wv_diff_ten[1],
			      double& mom_u_source,
			      double& dmom_u_source_p,
			      double& mom_v_source,
			      double& dmom_v_source_p,
			      double& mom_w_source,
			      double& dmom_w_source_p,
			      double& mom_u_ham,
			      double dmom_u_ham_grad_p[nSpace],
			      double& mom_v_ham,
			      double dmom_v_ham_grad_p[nSpace],
			      double& mom_w_ham,
			      double dmom_w_ham_grad_p[nSpace])
    {
      double rho,nu,mu,H_rho,d_rho,H_mu,d_mu,norm_n;
      H_rho = smoothedHeaviside(eps_rho,phi);
      d_rho = smoothedDirac(eps_rho,phi);
      H_mu = smoothedHeaviside(eps_mu,phi);
      d_mu = smoothedDirac(eps_mu,phi);
  
      rho = rho_0*(1.0-H_rho)+rho_1*H_rho;
      nu  = nu_0*(1.0-H_mu)+nu_1*H_mu;
      mu  = rho_0*nu_0*(1.0-H_mu)+rho_1*nu_1*H_mu;

      mass_acc=1.0-C;

      //u momentum accumulation
      mom_u_acc=u*(1.0-C);
      dmom_u_acc_u=(1.0-C);
  
      //v momentum accumulation
      mom_v_acc=v*(1.0-C);
      dmom_v_acc_v=(1.0-C);
  
      //w momentum accumulation
      mom_w_acc=w*(1.0-C);
      dmom_w_acc_w=(1.0-C);
  
  
      //mass advective flux
      mass_adv[0]=u*(1.0-C);
      mass_adv[1]=v*(1.0-C);
      mass_adv[2]=w*(1.0-C);
  
      dmass_adv_u[0]=(1.0-C);
      dmass_adv_u[1]=0.0;
      dmass_adv_u[2]=0.0;

      dmass_adv_v[0]=0.0;
      dmass_adv_v[1]=(1.0-C);
      dmass_adv_v[2]=0.0;

      dmass_adv_w[0]=0.0;
      dmass_adv_w[1]=0.0;
      dmass_adv_w[2]=(1.0-C);

      //u momentum advective flux
      mom_u_adv[0]=u*u*(1.0-C);
      mom_u_adv[1]=u*v*(1.0-C);
      mom_u_adv[2]=u*w*(1.0-C);
  
      dmom_u_adv_u[0]=2.0*u*(1.0-C);
      dmom_u_adv_u[1]=v*(1.0-C);
      dmom_u_adv_u[2]=w*(1.0-C);
  
      dmom_u_adv_v[0]=0.0;
      dmom_u_adv_v[1]=u*(1.0-C);
      dmom_u_adv_v[2]=0.0;
  
      dmom_u_adv_w[0]=0.0;
      dmom_u_adv_w[1]=0.0;
      dmom_u_adv_w[2]=u*(1.0-C);
  
      //v momentum advective_flux
      mom_v_adv[0]=v*u*(1.0-C);
      mom_v_adv[1]=v*v*(1.0-C);
      mom_v_adv[2]=v*w*(1.0-C);
  
      dmom_v_adv_u[0]=v*(1.0-C);
      dmom_v_adv_u[1]=0.0;
      dmom_v_adv_u[2]=0.0;
  
      dmom_v_adv_w[0]=0.0;
      dmom_v_adv_w[1]=0.0;
      dmom_v_adv_w[2]=v*(1.0-C);
  
      dmom_v_adv_v[0]=u*(1.0-C);
      dmom_v_adv_v[1]=2.0*v*(1.0-C);
      dmom_v_adv_v[2]=w*(1.0-C);
  
      //w momentum advective_flux
      mom_w_adv[0]=w*u*(1.0-C);
      mom_w_adv[1]=w*v*(1.0-C);
      mom_w_adv[2]=w*w*(1.0-C);
  
      dmom_w_adv_u[0]=w*(1.0-C);
      dmom_w_adv_u[1]=0.0;
      dmom_w_adv_u[2]=0.0;
  
      dmom_w_adv_v[0]=0.0;
      dmom_w_adv_v[1]=w*(1.0-C);
      dmom_w_adv_v[2]=0.0;
  
      dmom_w_adv_w[0]=u*(1.0-C);
      dmom_w_adv_w[1]=v*(1.0-C);
      dmom_w_adv_w[2]=2.0*w*(1.0-C);

      //u momentum diffusion tensor
      mom_u_diff_ten[0] = 2.0*nu*(1.0-C);
      mom_u_diff_ten[1] = nu*(1.0-C);
      mom_u_diff_ten[2] = nu*(1.0-C);
  
      mom_uv_diff_ten[0]=nu*(1.0-C);
  
      mom_uw_diff_ten[0]=nu*(1.0-C);
  
      //v momentum diffusion tensor
      mom_v_diff_ten[0] = nu*(1.0-C);
      mom_v_diff_ten[1] = 2.0*nu*(1.0-C);
      mom_v_diff_ten[2] = nu*(1.0-C);
  
      mom_vu_diff_ten[0]=nu*(1.0-C);
  
      mom_vw_diff_ten[0]=nu*(1.0-C);
  
      //w momentum diffusion tensor
      mom_w_diff_ten[0] = nu*(1.0-C);
      mom_w_diff_ten[1] = nu*(1.0-C);
      mom_w_diff_ten[2] = 2.0*nu*(1.0-C);
  
      mom_wu_diff_ten[0]=nu*(1.0-C);
  
      mom_wv_diff_ten[0]=nu*(1.0-C);
  
      //momentum sources
      norm_n = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      mom_u_source = -g[0]*(1.0-C) + sedimentForce[0] - p*grad_C[0]/rho;// - d_mu*sigma*kappa*n[0]/(rho*(norm_n+1.0e-8));
      dmom_u_source_p = -grad_C[0]/rho;
      mom_v_source = -g[1]*(1.0-C) + sedimentForce[1] - p*grad_C[1]/rho;// - d_mu*sigma*kappa*n[1]/(rho*(norm_n+1.0e-8));
      dmom_v_source_p = -grad_C[1]/rho;
      mom_w_source = -g[2]*(1.0-C) + sedimentForce[2] - p*grad_C[2]/rho;// - d_mu*sigma*kappa*n[2]/(rho*(norm_n+1.0e-8));
      dmom_w_source_p = -grad_C[2]/rho;
   
      //u momentum Hamiltonian (pressure)
      mom_u_ham = grad_p[0]*(1.0-C)/rho;
      dmom_u_ham_grad_p[0]=1.0*(1.0-C)/rho;
      dmom_u_ham_grad_p[1]=0.0;
      dmom_u_ham_grad_p[2]=0.0;
  
      //v momentum Hamiltonian (pressure)
      mom_v_ham = grad_p[1]*(1.0-C)/rho;
      dmom_v_ham_grad_p[0]=0.0;
      dmom_v_ham_grad_p[1]=1.0*(1.0-C)/rho;
      dmom_v_ham_grad_p[2]=0.0;
  
      //w momentum Hamiltonian (pressure)
      mom_w_ham = grad_p[2]*(1.0-C)/rho;
      dmom_w_ham_grad_p[0]=0.0;
      dmom_w_ham_grad_p[1]=0.0;
      dmom_w_ham_grad_p[2]=1.0*(1.0-C)/rho;
    }
  
    inline
    void calculateSubgridError_tau(const double&  hFactor,
				   const double& elementDiameter,
				   const double& dmt,
				   const double& dm,
				   const double df[nSpace],
				   const double& a,
				   const double&  pfac,
				   double& tau_v,
				   double& tau_p,
				   double& cfl)
    {
      double h,oneByAbsdt,density,viscosity,nrm_df;
      h = hFactor*elementDiameter;
      density = dm;
      viscosity =  a;
      nrm_df=0.0;
      for(int I=0;I<nSpace;I++)
	nrm_df+=df[I]*df[I];
      nrm_df = sqrt(nrm_df);
      cfl = nrm_df/(h*density);//this is really cfl/dt, but that's what we want to know, the step controller expect this
      oneByAbsdt =  fabs(dmt);
      tau_v = 1.0/(4.0*viscosity/(h*h) + 2.0*nrm_df/h + oneByAbsdt);
      tau_p = (4.0*viscosity + 2.0*nrm_df*h + oneByAbsdt*h*h)/pfac;
    }


    inline
    void calculateSubgridError_tau(     const double&  Ct_sge,
                                        const double&  Cd_sge,
			                const double   G[nSpace*nSpace],
					const double&  G_dd_G,
					const double&  tr_G,
					const double&  A0,
					const double   Ai[nSpace],
					const double&  Kij,
					const double&  pfac,
					double& tau_v,
					double& tau_p,
					double& q_cfl)	
    {
      double v_d_Gv=0.0; 
      for(int I=0;I<nSpace;I++) 
         for (int J=0;J<nSpace;J++) 
           v_d_Gv += Ai[I]*G[I*nSpace+J]*Ai[J];     
    
      tau_v = 1.0/sqrt(Ct_sge*A0*A0 + v_d_Gv + Cd_sge*Kij*Kij*G_dd_G); 
      tau_p = 1.0/(pfac*tr_G*tau_v);     
    }



    /* inline */
    /*   void calculateSubgridError_tau(const double&  Ct_sge, */
    /* 				   const double&  Cd_sge, */
    /* 				   const double* G, */
    /* 				   const double& G_dd_G, */
    /* 				   const double& tr_G, */
    /* 				   const double& rho, */
    /* 				   const double& Dt, */
    /* 				   const double v[nSpace], */
    /* 				   const double& mu, */
    /* 				   double& tau_p, */
    /* 				   double& tau_v, */
    /* 				   double& cfl) */
    /* { */
    /*   const double rho2=rho*rho,Dt2=Dt*Dt,mu2=mu*mu; */
    /*   register double v_d_Gv=0.0; */
    /*   for(int I=0;I<nSpace;I++) */
    /*     for (int J=0;J<nSpace;J++) */
    /* 	v_d_Gv += v[I]*G[I*nSpace+J]*v[J]; */
    /*   cfl = 2.0*sqrt(v_d_Gv); */
    /*   //cek 1.0/sqrt(rho2*Dt2 + 4*v_d_Gv + 144*mu2*G_dd_G); ? */
    /*   /\* tau_v = 1.0/sqrt(rho2*Dt2 + 4*v_d_Gv + 144*mu2*G_dd_G); *\/ */
    /*   /\* tau_p = 1.0/(tr_G*tau_v);  *\/ */
    /*   //cek "correct" tau */
    /*   tau_v = 1.0/sqrt(Ct_sge*rho2*Dt2 + rho2*v_d_Gv + Cd_sge*mu2*G_dd_G); */
    /*   tau_p = 1.0/(tr_G*tau_v); */
    /*   //debug */
    /*   /\* double tau_v_old = tau_v,tau_p_old = tau_p; *\/ */
    /*   /\* double nrm_v=0.0,h=1.0/20.0; *\/ */
    /*   /\* double oneByAbsdt =  fabs(Dt); *\/ */
    /*   /\* for(int I=0;I<nSpace;I++) *\/ */
    /*   /\* 	nrm_v += v[I]*v[I]; *\/ */
    /*   /\* nrm_v = sqrt(nrm_v); *\/ */
    /*   /\* cfl = nrm_v/h; *\/ */
    /*   /\* tau_v = 1.0/(4.0*mu/(h*h) + 2.0*rho*nrm_v/h + oneByAbsdt); *\/ */
    /*   /\* tau_p = 4.0*mu + 2.0*rho*nrm_v*h+ oneByAbsdt*h*h; *\/ */
    /*   /\* std::cout<<"nrm_v "<<nrm_v<<" tau_v "<<tau_v<<"\t"<<tau_v_old<<" tau_p "<<tau_p<<'\t'<<tau_p_old<<std::endl; *\/ */
    /* } */

    inline
    void calculateSubgridError_tauRes(const double& tau_p,
				      const double& tau_v,
				      const double& pdeResidualP,
				      const double& pdeResidualU,
				      const double& pdeResidualV,
				      const double& pdeResidualW,
				      double& subgridErrorP,
				      double& subgridErrorU,
				      double& subgridErrorV,
				      double& subgridErrorW)
    {
      /* GLS pressure */
      subgridErrorP = -tau_p*pdeResidualP;
      /* GLS momentum */
      subgridErrorU = -tau_v*pdeResidualU;
      subgridErrorV = -tau_v*pdeResidualV;
      subgridErrorW = -tau_v*pdeResidualW;
    }

    inline
    void calculateSubgridErrorDerivatives_tauRes(const double& tau_p,
						 const double& tau_v,
						 const double dpdeResidualP_du[nDOF_trial_element],
						 const double dpdeResidualP_dv[nDOF_trial_element],
						 const double dpdeResidualP_dw[nDOF_trial_element],
						 const double dpdeResidualU_dp[nDOF_trial_element],
						 const double dpdeResidualU_du[nDOF_trial_element],
						 const double dpdeResidualV_dp[nDOF_trial_element],
						 const double dpdeResidualV_dv[nDOF_trial_element],
						 const double dpdeResidualW_dp[nDOF_trial_element],
						 const double dpdeResidualW_dw[nDOF_trial_element],
						 double dsubgridErrorP_du[nDOF_trial_element],
						 double dsubgridErrorP_dv[nDOF_trial_element],
						 double dsubgridErrorP_dw[nDOF_trial_element],
						 double dsubgridErrorU_dp[nDOF_trial_element],
						 double dsubgridErrorU_du[nDOF_trial_element],
						 double dsubgridErrorV_dp[nDOF_trial_element],
						 double dsubgridErrorV_dv[nDOF_trial_element],
						 double dsubgridErrorW_dp[nDOF_trial_element],
						 double dsubgridErrorW_dw[nDOF_trial_element])
    {
      for (int j=0;j<nDOF_trial_element;j++)
	{
	  /* GLS pressure */
	  dsubgridErrorP_du[j] = -tau_p*dpdeResidualP_du[j];
	  dsubgridErrorP_dv[j] = -tau_p*dpdeResidualP_dv[j];
	  dsubgridErrorP_dw[j] = -tau_p*dpdeResidualP_dw[j];
	  /* GLS  momentum*/
	  /* u */
	  dsubgridErrorU_dp[j] = -tau_v*dpdeResidualU_dp[j];
	  dsubgridErrorU_du[j] = -tau_v*dpdeResidualU_du[j];
	  /* v */
	  dsubgridErrorV_dp[j] = -tau_v*dpdeResidualV_dp[j];
	  dsubgridErrorV_dv[j] = -tau_v*dpdeResidualV_dv[j];
	  /* w */
	  dsubgridErrorW_dp[j] = -tau_v*dpdeResidualW_dp[j];
	  dsubgridErrorW_dw[j] = -tau_v*dpdeResidualW_dw[j];
	}
    }

    inline
    void exteriorNumericalAdvectiveFlux(const int& isDOFBoundary_p,
					const int& isDOFBoundary_u,
					const int& isDOFBoundary_v,
					const int& isDOFBoundary_w,
					const int& isFluxBoundary_p,
					const int& isFluxBoundary_u,
					const int& isFluxBoundary_v,
					const int& isFluxBoundary_w,
					const double& oneByRho,
					const double n[nSpace],
					const double& bc_p,
					const double bc_f_mass[nSpace],
					const double bc_f_umom[nSpace],
					const double bc_f_vmom[nSpace],
					const double bc_f_wmom[nSpace],
					const double& bc_flux_mass,
					const double& bc_flux_umom,
					const double& bc_flux_vmom,
					const double& bc_flux_wmom,
					const double& p,
					const double f_mass[nSpace],
					const double f_umom[nSpace],
					const double f_vmom[nSpace],
					const double f_wmom[nSpace],
					const double df_mass_du[nSpace],
					const double df_mass_dv[nSpace],
					const double df_mass_dw[nSpace],
					const double df_umom_dp[nSpace],
					const double df_umom_du[nSpace],
					const double df_umom_dv[nSpace],
					const double df_umom_dw[nSpace],
					const double df_vmom_dp[nSpace],
					const double df_vmom_du[nSpace],
					const double df_vmom_dv[nSpace],
					const double df_vmom_dw[nSpace],
					const double df_wmom_dp[nSpace],
					const double df_wmom_du[nSpace],
					const double df_wmom_dv[nSpace],
					const double df_wmom_dw[nSpace],
					double& flux_mass,
					double& flux_umom,
					double& flux_vmom,
					double& flux_wmom,
					double* velocity)
    {
      double flowDirection;
      flux_mass = 0.0;
      flux_umom = 0.0;
      flux_vmom = 0.0;
      flux_wmom = 0.0;
      flowDirection=n[0]*f_mass[0]+n[1]*f_mass[1]+n[2]*f_mass[2];
      if (isDOFBoundary_u != 1)
	{
	  flux_mass += n[0]*f_mass[0];
	  velocity[0] = f_mass[0];
	  if (flowDirection >= 0.0)
	    {
	      flux_umom += n[0]*f_umom[0];
	      flux_vmom += n[0]*f_vmom[0];
	      flux_wmom += n[0]*f_wmom[0];
	    }
	}
      else
	{
	  flux_mass += n[0]*bc_f_mass[0];
	  velocity[0] = bc_f_mass[0];
	  //cek still upwind the advection for Dirichlet?
	  if (flowDirection >= 0.0)
	    {
	      flux_umom += n[0]*f_umom[0];
	      flux_vmom += n[0]*f_vmom[0];
	      flux_wmom += n[0]*f_wmom[0];
	    }
	  else
	    {
	      flux_umom+=n[0]*bc_f_umom[0];
	      flux_vmom+=n[0]*bc_f_vmom[0];
	      flux_wmom+=n[0]*bc_f_wmom[0];
	    }
	}
      if (isDOFBoundary_v != 1)
	{
	  flux_mass+=n[1]*f_mass[1];
	  velocity[1] = f_mass[1];
	  if (flowDirection >= 0.0)
	    {
	      flux_umom+=n[1]*f_umom[1];
	      flux_vmom+=n[1]*f_vmom[1];
	      flux_wmom+=n[1]*f_wmom[1];
	    }
	}
      else
	{
	  flux_mass+=n[1]*bc_f_mass[1];
	  velocity[1] = bc_f_mass[1];
	  //cek still upwind the advection for Dirichlet?
	  if (flowDirection >= 0.0)
	    {
	      flux_umom+=n[1]*f_umom[1];
	      flux_vmom+=n[1]*f_vmom[1];
	      flux_wmom+=n[1]*f_wmom[1];
	    }
	  else
	    {
	      flux_umom+=n[1]*bc_f_umom[1];
	      flux_vmom+=n[1]*bc_f_vmom[1];
	      flux_wmom+=n[1]*bc_f_wmom[1];
	    }
	}
      if (isDOFBoundary_w != 1)
	{
	  flux_mass+=n[2]*f_mass[2];
	  velocity[2] = f_mass[2];
	  if (flowDirection >= 0.0)
	    {
	      flux_umom+=n[2]*f_umom[2];
	      flux_vmom+=n[2]*f_vmom[2];
	      flux_wmom+=n[2]*f_wmom[2];
	    }
	}
      else
	{
	  flux_mass +=n[2]*bc_f_mass[2];
	  velocity[2] = bc_f_mass[2];
	  //cek still upwind the advection for Dirichlet?
	  if (flowDirection >= 0.0)
	    {
	      flux_umom+=n[2]*f_umom[2];
	      flux_vmom+=n[2]*f_vmom[2];
	      flux_wmom+=n[2]*f_wmom[2];
	    }
	  else
	    {
	      flux_umom+=n[2]*bc_f_umom[2];
	      flux_vmom+=n[2]*bc_f_vmom[2];
	      flux_wmom+=n[2]*bc_f_wmom[2];
	    }
	}
      if (isDOFBoundary_p == 1)
	{
	  flux_umom+= n[0]*(bc_p-p)*oneByRho;
	  flux_vmom+= n[1]*(bc_p-p)*oneByRho;
	  flux_wmom+= n[2]*(bc_p-p)*oneByRho;
	}
      if (isFluxBoundary_p == 1)
	{
	  velocity[0] += (bc_flux_mass - flux_mass)*n[0];
	  velocity[1] += (bc_flux_mass - flux_mass)*n[1];
	  velocity[2] += (bc_flux_mass - flux_mass)*n[2];
	  flux_mass = bc_flux_mass;
	}
      if (isFluxBoundary_u == 1)
	{
	  flux_umom = bc_flux_umom;
	}
      if (isFluxBoundary_v == 1)
	{
	  flux_vmom = bc_flux_vmom;
	}
      if (isFluxBoundary_w == 1)
	{
	  flux_wmom = bc_flux_wmom;
	}
    }

    inline
    void exteriorNumericalAdvectiveFluxDerivatives(const int& isDOFBoundary_p,
						   const int& isDOFBoundary_u,
						   const int& isDOFBoundary_v,
						   const int& isDOFBoundary_w,
						   const int& isFluxBoundary_p,
						   const int& isFluxBoundary_u,
						   const int& isFluxBoundary_v,
						   const int& isFluxBoundary_w,
						   const double& oneByRho,
						   const double n[nSpace],
						   const double& bc_p,
						   const double bc_f_mass[nSpace],
						   const double bc_f_umom[nSpace],
						   const double bc_f_vmom[nSpace],
						   const double bc_f_wmom[nSpace],
						   const double& bc_flux_mass,
						   const double& bc_flux_umom,
						   const double& bc_flux_vmom,
						   const double& bc_flux_wmom,
						   const double& p,
						   const double f_mass[nSpace],
						   const double f_umom[nSpace],
						   const double f_vmom[nSpace],
						   const double f_wmom[nSpace],
						   const double df_mass_du[nSpace],
						   const double df_mass_dv[nSpace],
						   const double df_mass_dw[nSpace],
						   const double df_umom_dp[nSpace],
						   const double df_umom_du[nSpace],
						   const double df_umom_dv[nSpace],
						   const double df_umom_dw[nSpace],
						   const double df_vmom_dp[nSpace],
						   const double df_vmom_du[nSpace],
						   const double df_vmom_dv[nSpace],
						   const double df_vmom_dw[nSpace],
						   const double df_wmom_dp[nSpace],
						   const double df_wmom_du[nSpace],
						   const double df_wmom_dv[nSpace],
						   const double df_wmom_dw[nSpace],
						   double& dflux_mass_du,
						   double& dflux_mass_dv,
						   double& dflux_mass_dw,
						   double& dflux_umom_dp,
						   double& dflux_umom_du,
						   double& dflux_umom_dv,
						   double& dflux_umom_dw,
						   double& dflux_vmom_dp,
						   double& dflux_vmom_du,
						   double& dflux_vmom_dv,
						   double& dflux_vmom_dw,
						   double& dflux_wmom_dp,
						   double& dflux_wmom_du,
						   double& dflux_wmom_dv,
						   double& dflux_wmom_dw)
    {
      double flowDirection;
      dflux_mass_du = 0.0;
      dflux_mass_dv = 0.0;
      dflux_mass_dw = 0.0;
  
      dflux_umom_dp = 0.0;
      dflux_umom_du = 0.0;
      dflux_umom_dv = 0.0;
      dflux_umom_dw = 0.0;
  
      dflux_vmom_dp = 0.0;
      dflux_vmom_du = 0.0;
      dflux_vmom_dv = 0.0;
      dflux_vmom_dw = 0.0;
  
      dflux_wmom_dp = 0.0;
      dflux_wmom_du = 0.0;
      dflux_wmom_dv = 0.0;
      dflux_wmom_dw = 0.0;
  
      flowDirection=n[0]*f_mass[0]+n[1]*f_mass[1]+n[2]*f_mass[2];
      if (isDOFBoundary_u != 1)
	{
	  dflux_mass_du += n[0]*df_mass_du[0];
	  if (flowDirection >= 0.0)
	    {
	      dflux_umom_du += n[0]*df_umom_du[0];
	      dflux_vmom_du += n[0]*df_vmom_du[0];
	      dflux_vmom_dv += n[0]*df_vmom_dv[0];
	      dflux_wmom_du += n[0]*df_wmom_du[0];
	      dflux_wmom_dw += n[0]*df_wmom_dw[0];
	    }
	}
      else
	{
	  //cek still upwind the advection for Dirichlet?
	  if (flowDirection >= 0.0)
	    {
	      dflux_umom_du += n[0]*df_umom_du[0];
	      dflux_vmom_du += n[0]*df_vmom_du[0];
	      dflux_vmom_dv += n[0]*df_vmom_dv[0];
	      dflux_wmom_du += n[0]*df_wmom_du[0];
	      dflux_wmom_dw += n[0]*df_wmom_dw[0];
	    }
	  else
	    {
	      if (isDOFBoundary_v != 1)
		dflux_vmom_dv += n[0]*df_vmom_dv[0];
	      if (isDOFBoundary_w != 1)
		dflux_wmom_dw += n[0]*df_wmom_dw[0];
	    }
	}
      if (isDOFBoundary_v != 1)
	{
	  dflux_mass_dv += n[1]*df_mass_dv[1];
	  if (flowDirection >= 0.0)
	    {
	      dflux_umom_du += n[1]*df_umom_du[1];
	      dflux_umom_dv += n[1]*df_umom_dv[1];
	      dflux_vmom_dv += n[1]*df_vmom_dv[1];
	      dflux_wmom_dw += n[1]*df_wmom_dw[1];
	      dflux_wmom_dv += n[1]*df_wmom_dv[1];
	    }
	}
      else
	{
	  //cek still upwind the advection for Dirichlet?
	  if (flowDirection >= 0.0)
	    {
	      dflux_umom_du += n[1]*df_umom_du[1];
	      dflux_umom_dv += n[1]*df_umom_dv[1];
	      dflux_vmom_dv += n[1]*df_vmom_dv[1];
	      dflux_wmom_dw += n[1]*df_wmom_dw[1];
	      dflux_wmom_dv += n[1]*df_wmom_dv[1];
	    }
	  else
	    {
	      if (isDOFBoundary_u != 1)
		dflux_umom_du += n[1]*df_umom_du[1];
	      if (isDOFBoundary_w != 1)
		dflux_wmom_dw += n[1]*df_wmom_dw[1];
	    }
	}
      if (isDOFBoundary_w != 1)
	{
	  dflux_mass_dw+=n[2]*df_mass_dw[2];
	  if (flowDirection >= 0.0)
	    {
	      dflux_umom_du += n[2]*df_umom_du[2];
	      dflux_umom_dw += n[2]*df_umom_dw[2];
	      dflux_vmom_dv += n[2]*df_vmom_dv[2];
	      dflux_vmom_dw += n[2]*df_vmom_dw[2];
	      dflux_wmom_dw += n[2]*df_wmom_dw[2];
	    }
	}
      else
	{
	  //cek still upwind the advection for Dirichlet?
	  if (flowDirection >= 0.0)
	    {
	      dflux_umom_du += n[2]*df_umom_du[2];
	      dflux_umom_dw += n[2]*df_umom_dw[2];
	      dflux_vmom_dv += n[2]*df_vmom_dv[2];
	      dflux_vmom_dw += n[2]*df_vmom_dw[2];
	      dflux_wmom_dw += n[2]*df_wmom_dw[2];
	    }
	  else
	    {
	      if (isDOFBoundary_u != 1)
		dflux_umom_du += n[2]*df_umom_du[2];
	      if (isDOFBoundary_v != 1)
		dflux_vmom_dv += n[2]*df_vmom_dv[2];
	    }
	}
      if (isDOFBoundary_p == 1)
	{
	  dflux_umom_dp= -n[0]*oneByRho;
	  dflux_vmom_dp= -n[1]*oneByRho;
	  dflux_wmom_dp= -n[2]*oneByRho;
	}
      if (isFluxBoundary_p == 1)
	{
	  dflux_mass_du = 0.0;
	  dflux_mass_dv = 0.0;
	  dflux_mass_dw = 0.0;
	}
      if (isFluxBoundary_u == 1)
	{
	  dflux_umom_dp = 0.0;
	  dflux_umom_du = 0.0;
	  dflux_umom_dv = 0.0;
	  dflux_umom_dw = 0.0;
	}
      if (isFluxBoundary_v == 1)
	{
	  dflux_vmom_dp = 0.0;
	  dflux_vmom_du = 0.0;
	  dflux_vmom_dv = 0.0;
	  dflux_vmom_dw = 0.0;
	}
      if (isFluxBoundary_w == 1)
	{
	  dflux_wmom_dp = 0.0;
	  dflux_wmom_du = 0.0;
	  dflux_wmom_dv = 0.0;
	  dflux_wmom_dw = 0.0;
	}
    }

    inline
    void exteriorNumericalDiffusiveFlux(const double& eps,
					const double& phi,
					int* rowptr,
					int* colind,
					const int& isDOFBoundary,
					const int& isFluxBoundary,
					const double n[nSpace],
					double* bc_a,
					const double& bc_u,
					const double& bc_flux,
					double* a,
					const double grad_phi[nSpace],
					const double& u,
					const double& penalty,
					double& flux)
    {
      double diffusiveVelocityComponent_I,penaltyFlux,max_a;
      if(isDOFBoundary == 1)
	{
	  flux = 0.0;
	  max_a=0.0;
	  for(int I=0;I<nSpace;I++)
	    {
	      diffusiveVelocityComponent_I=0.0;
	      for(int m=rowptr[I];m<rowptr[I+1];m++)
		{
		  diffusiveVelocityComponent_I -= a[m]*grad_phi[colind[m]];
		  max_a = fmax(max_a,a[m]);
		}
	      flux+= diffusiveVelocityComponent_I*n[I];
	    }
	  penaltyFlux = max_a*penalty*(u-bc_u);
	  flux += penaltyFlux;
	  //contact line slip
	  //flux*=(smoothedDirac(eps,0) - smoothedDirac(eps,phi))/smoothedDirac(eps,0);
	}
      else if(isFluxBoundary == 1)
	{
	  flux = bc_flux;
	}
      else
	{
	  std::cerr<<"warning, diffusion term with no boundary condition set, setting diffusive flux to 0.0"<<std::endl;
	  flux = 0.0;
	}
    }


    inline
    double ExteriorNumericalDiffusiveFluxJacobian(const double& eps,
						  const double& phi,
						  int* rowptr,
						  int* colind,
						  const int& isDOFBoundary,
						  const double n[nSpace],
						  double* a,
						  const double& v,
						  const double grad_v[nSpace],
						  const double& penalty)
    {
      double dvel_I,tmp=0.0,max_a=0.0;
      if(isDOFBoundary >= 1)
	{
	  for(int I=0;I<nSpace;I++)
	    {
	      dvel_I=0.0;
	      for(int m=rowptr[I];m<rowptr[I+1];m++)
		{
		  dvel_I -= a[m]*grad_v[colind[m]];
		  max_a = fmax(max_a,a[m]);
		}
	      tmp += dvel_I*n[I];
	    }
	  tmp +=max_a*penalty*v;
	  //contact line slip
	  //tmp*=(smoothedDirac(eps,0) - smoothedDirac(eps,phi))/smoothedDirac(eps,0);
	}
      return tmp;
    }

    void calculateResidual(//element
			   double* mesh_trial_ref,
			   double* mesh_grad_trial_ref,
			   double* mesh_dof,
			   double* mesh_velocity_dof,
			   double MOVING_DOMAIN,
			   int* mesh_l2g,
			   double* dV_ref,
			   double* p_trial_ref,
			   double* p_grad_trial_ref,
			   double* p_test_ref,
			   double* p_grad_test_ref,
			   double* vel_trial_ref,
			   double* vel_grad_trial_ref,
			   double* vel_test_ref,
			   double* vel_grad_test_ref,
			   //element boundary
			   double* mesh_trial_trace_ref,
			   double* mesh_grad_trial_trace_ref,
			   double* dS_ref,
			   double* p_trial_trace_ref,
			   double* p_grad_trial_trace_ref,
			   double* p_test_trace_ref,
			   double* p_grad_test_trace_ref,
			   double* vel_trial_trace_ref,
			   double* vel_grad_trial_trace_ref,
			   double* vel_test_trace_ref,
			   double* vel_grad_test_trace_ref,					 
			   double* normal_ref,
			   double* boundaryJac_ref,
			   //physics
			   double* elementDiameter,
			   double hFactor,
			   int nElements_global,
			   double useRBLES,
			   double useMetrics, 
			   double alphaBDF,
			   double epsFact_rho,
			   double epsFact_mu, 
			   double sigma,
			   double rho_0,
			   double nu_0,
			   double rho_1,
			   double nu_1,
			   double Ct_sge,
			   double Cd_sge,
			   double C_dc,
			   int* p_l2g, 
			   int* vel_l2g, 
			   double* C_dof, 
			   double* p_dof, 
			   double* u_dof, 
			   double* v_dof, 
			   double* w_dof,
			   double* g,
			   double* q_sedimentForce,
			   //double* ebqe_sedimentForce,
			   double* phi,
			   double* normal_phi,
			   double* kappa_phi,
			   double* q_mass_acc,
			   double* q_mom_u_acc,
			   double* q_mom_v_acc,
			   double* q_mom_w_acc,
			   double* q_mass_adv,
			   double* q_mass_acc_beta_bdf,double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
			   double* q_velocity_sge,
			   double* q_cfl,
			   double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
			   double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
			   int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
			   int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
			   int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
			   int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
			   int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
			   int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
			   int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
			   int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
			   int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
			   int offset_p, int offset_u, int offset_v, int offset_w, 
			   int stride_p, int stride_u, int stride_v, int stride_w, 
			   double* globalResidual,
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
			   double* ebqe_phi_ext,
			   double* ebqe_normal_phi_ext,
			   double* ebqe_kappa_phi_ext,
			   int* isDOFBoundary_p,
			   int* isDOFBoundary_u,
			   int* isDOFBoundary_v,
			   int* isDOFBoundary_w,
			   int* isAdvectiveFluxBoundary_p,
			   int* isAdvectiveFluxBoundary_u,
			   int* isAdvectiveFluxBoundary_v,
			   int* isAdvectiveFluxBoundary_w,
			   int* isDiffusiveFluxBoundary_u,
			   int* isDiffusiveFluxBoundary_v,
			   int* isDiffusiveFluxBoundary_w,
			   double* ebqe_bc_p_ext,
			   double* ebqe_bc_flux_mass_ext,
			   double* ebqe_bc_flux_mom_u_adv_ext,
			   double* ebqe_bc_flux_mom_v_adv_ext,
			   double* ebqe_bc_flux_mom_w_adv_ext,
			   double* ebqe_bc_u_ext,
			   double* ebqe_bc_flux_u_diff_ext,
			   double* ebqe_penalty_ext,
			   double* ebqe_bc_v_ext,
			   double* ebqe_bc_flux_v_diff_ext,
			   double* ebqe_bc_w_ext,
			   double* ebqe_bc_flux_w_diff_ext,
			   double* q_velocity,
			   double* ebqe_velocity,
			   double* flux,
			   double* elementResidual_p_save)
    {
      //
      //loop over elements to compute volume integrals and load them into element and global residual
      //
      double globalConservationError=0.0;
      for(int eN=0;eN<nElements_global;eN++)
	{
	  //declare local storage for element residual and initialize
	  register double elementResidual_p[nDOF_test_element],
	    elementResidual_u[nDOF_test_element],
	    elementResidual_v[nDOF_test_element],
	    elementResidual_w[nDOF_test_element],
	    eps_rho,eps_mu;
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      elementResidual_p_save[eN_i]=0.0;
	      elementResidual_p[i]=0.0;
	      elementResidual_u[i]=0.0;
	      elementResidual_v[i]=0.0;
	      elementResidual_w[i]=0.0;
	    }//i
	  //
	  //loop over quadrature points and compute integrands
	  //
	  for(int k=0;k<nQuadraturePoints_element;k++)
	    {
	      //compute indices and declare local storage
	      register int eN_k = eN*nQuadraturePoints_element+k,
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element;
	      register double p=0.0,C=0.0,u=0.0,v=0.0,w=0.0,
		grad_p[nSpace],grad_C[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
		mass_acc=0.0,
		dmass_acc=0.0,
		mom_u_acc=0.0,
		dmom_u_acc_u=0.0,
		mom_v_acc=0.0,
		dmom_v_acc_v=0.0,
		mom_w_acc=0.0,
		dmom_w_acc_w=0.0,
		mass_adv[nSpace],
		dmass_adv_u[nSpace],
		dmass_adv_v[nSpace],
		dmass_adv_w[nSpace],
		mom_u_adv[nSpace],
		dmom_u_adv_u[nSpace],
		dmom_u_adv_v[nSpace],
		dmom_u_adv_w[nSpace],
		mom_v_adv[nSpace],
		dmom_v_adv_u[nSpace],
		dmom_v_adv_v[nSpace],
		dmom_v_adv_w[nSpace],
		mom_w_adv[nSpace],
		dmom_w_adv_u[nSpace],
		dmom_w_adv_v[nSpace],
		dmom_w_adv_w[nSpace],
		mom_u_diff_ten[nSpace],
		mom_v_diff_ten[nSpace],
		mom_w_diff_ten[nSpace],
		mom_uv_diff_ten[1],
		mom_uw_diff_ten[1],
		mom_vu_diff_ten[1],
		mom_vw_diff_ten[1],
		mom_wu_diff_ten[1],
		mom_wv_diff_ten[1],
		mom_u_source=0.0,
		dmom_u_source_p=0.0,
		mom_v_source=0.0,
		dmom_v_source_p=0.0,
		mom_w_source=0.0,
		dmom_w_source_p=0.0,
		mom_u_ham=0.0,
		dmom_u_ham_grad_p[nSpace],
		mom_v_ham=0.0,
		dmom_v_ham_grad_p[nSpace],
		mom_w_ham=0.0,
		dmom_w_ham_grad_p[nSpace],
		mass_acc_t=0.0,
		dmass_acc_t=0.0,
		mom_u_acc_t=0.0,
		dmom_u_acc_u_t=0.0,
		mom_v_acc_t=0.0,
		dmom_v_acc_v_t=0.0,
		mom_w_acc_t=0.0,
		dmom_w_acc_w_t=0.0,
		pdeResidual_p=0.0,
		pdeResidual_u=0.0,
		pdeResidual_v=0.0,
		pdeResidual_w=0.0,
		Lstar_u_p[nDOF_test_element],
		Lstar_v_p[nDOF_test_element],
		Lstar_w_p[nDOF_test_element],
		Lstar_u_u[nDOF_test_element],
		Lstar_v_v[nDOF_test_element],
		Lstar_w_w[nDOF_test_element],
		Lstar_p_u[nDOF_test_element],
		Lstar_p_v[nDOF_test_element],
		Lstar_p_w[nDOF_test_element],
		subgridError_p=0.0,
		subgridError_u=0.0,
		subgridError_v=0.0,
		subgridError_w=0.0,
		tau_p=0.0,tau_p0=0.0,tau_p1=0.0,
		tau_v=0.0,tau_v0=0.0,tau_v1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
		p_test_dV[nDOF_trial_element],vel_test_dV[nDOF_trial_element],
		p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
		dV,x,y,z,xt,yt,zt,
		G[nSpace*nSpace],G_dd_G,tr_G,norm_Rv,h_phi, dmom_adv_star[nSpace],dmom_adv_sge[nSpace];
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
	      //xt=0.0;yt=0.0;zt=0.0;
	      //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<zt<<std::endl;
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.calculateG(jacInv,G,G_dd_G,tr_G);
	      ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);
	      
	      eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	     
	      //get the trial function gradients
	      ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      //get the solution
	      ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p);
	      ck.valFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],C);
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],u);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],v);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w);
	      //get the solution gradients
	      ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p);
	      ck.gradFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_C);
	      ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_u);
	      ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_v);
	      ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV;
		  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		      vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		    }
		}
	      //save velocity at quadrature points for other models to use
	      q_velocity[eN_k_nSpace+0]=u;
	      q_velocity[eN_k_nSpace+1]=v;
	      q_velocity[eN_k_nSpace+2]=w;
	      //
	      //calculate pde coefficients at quadrature points
	      //
	      //cek debug
	      //eps_rho=0.1;
	      //eps_mu=0.1;
	      evaluateCoefficients(eps_rho,
				   eps_mu,
				   sigma,
				   rho_0,
				   nu_0,
				   rho_1,
				   nu_1,
				   g,
				   &q_sedimentForce[eN_k_nSpace],
				   phi[eN_k],
				   &normal_phi[eN_k_nSpace],
				   kappa_phi[eN_k],
				   C,
				   grad_C,
				   p,
				   grad_p,
				   u,
				   v,
				   w,
				   mass_acc,
				   mom_u_acc,
				   dmom_u_acc_u,
				   mom_v_acc,
				   dmom_v_acc_v,
				   mom_w_acc,
				   dmom_w_acc_w,
				   mass_adv,
				   dmass_adv_u,
				   dmass_adv_v,
				   dmass_adv_w,
				   mom_u_adv,
				   dmom_u_adv_u,
				   dmom_u_adv_v,
				   dmom_u_adv_w,
				   mom_v_adv,
				   dmom_v_adv_u,
				   dmom_v_adv_v,
				   dmom_v_adv_w,
				   mom_w_adv,
				   dmom_w_adv_u,
				   dmom_w_adv_v,
				   dmom_w_adv_w,
				   mom_u_diff_ten,
				   mom_v_diff_ten,
				   mom_w_diff_ten,
				   mom_uv_diff_ten,
				   mom_uw_diff_ten,
				   mom_vu_diff_ten,
				   mom_vw_diff_ten,
				   mom_wu_diff_ten,
				   mom_wv_diff_ten,
				   mom_u_source,
				   dmom_u_source_p,
				   mom_v_source,
				   dmom_v_source_p,
				   mom_w_source,
				   dmom_w_source_p,
				   mom_u_ham,
				   dmom_u_ham_grad_p,
				   mom_v_ham,
				   dmom_v_ham_grad_p,
				   mom_w_ham,
				   dmom_w_ham_grad_p);          
	      //
	      //save momentum for time history and velocity for subgrid error
	      //
	      q_mass_acc[eN_k] = mass_acc;                            
	      q_mom_u_acc[eN_k] = mom_u_acc;                            
	      q_mom_v_acc[eN_k] = mom_v_acc;                            
	      q_mom_w_acc[eN_k] = mom_w_acc;
	      //subgrid error uses grid scale velocity
	      q_mass_adv[eN_k_nSpace+0] = u;
	      q_mass_adv[eN_k_nSpace+1] = v;
	      q_mass_adv[eN_k_nSpace+2] = w;
	      //
	      //moving mesh
	      //
	      //transform the continuity equation as if the accumulation term was  d(1)/dt
	      mass_adv[0] -= MOVING_DOMAIN*mass_acc*xt;
	      mass_adv[1] -= MOVING_DOMAIN*mass_acc*yt;
	      mass_adv[2] -= MOVING_DOMAIN*mass_acc*zt;

	      mom_u_adv[0] -= MOVING_DOMAIN*mom_u_acc*xt;
	      mom_u_adv[1] -= MOVING_DOMAIN*mom_u_acc*yt;
	      mom_u_adv[2] -= MOVING_DOMAIN*mom_u_acc*zt;
	      dmom_u_adv_u[0] -= MOVING_DOMAIN*dmom_u_acc_u*xt;
	      dmom_u_adv_u[1] -= MOVING_DOMAIN*dmom_u_acc_u*yt;
	      dmom_u_adv_u[2] -= MOVING_DOMAIN*dmom_u_acc_u*zt;

	      mom_v_adv[0] -= MOVING_DOMAIN*mom_v_acc*xt;
	      mom_v_adv[1] -= MOVING_DOMAIN*mom_v_acc*yt;
	      mom_v_adv[2] -= MOVING_DOMAIN*mom_v_acc*zt;
	      dmom_v_adv_v[0] -= MOVING_DOMAIN*dmom_v_acc_v*xt;
	      dmom_v_adv_v[1] -= MOVING_DOMAIN*dmom_v_acc_v*yt;
	      dmom_v_adv_v[2] -= MOVING_DOMAIN*dmom_v_acc_v*zt;

	      mom_w_adv[0] -= MOVING_DOMAIN*mom_w_acc*xt;
	      mom_w_adv[1] -= MOVING_DOMAIN*mom_w_acc*yt;
	      mom_w_adv[2] -= MOVING_DOMAIN*mom_w_acc*zt;
	      dmom_w_adv_w[0] -= MOVING_DOMAIN*dmom_w_acc_w*xt;
	      dmom_w_adv_w[1] -= MOVING_DOMAIN*dmom_w_acc_w*yt;
	      dmom_w_adv_w[2] -= MOVING_DOMAIN*dmom_w_acc_w*zt;
	      //
	      //calculate time derivative at quadrature points
	      //
	      ck.bdf(alphaBDF,
		     q_mass_acc_beta_bdf[eN_k],
		     mass_acc,
		     dmass_acc,
		     mass_acc_t,
		     dmass_acc_t);
	      ck.bdf(alphaBDF,
		     q_mom_u_acc_beta_bdf[eN_k],
		     mom_u_acc,
		     dmom_u_acc_u,
		     mom_u_acc_t,
		     dmom_u_acc_u_t);
	      ck.bdf(alphaBDF,
		     q_mom_v_acc_beta_bdf[eN_k],
		     mom_v_acc,
		     dmom_v_acc_v,
		     mom_v_acc_t,
		     dmom_v_acc_v_t);
	      ck.bdf(alphaBDF,
		     q_mom_w_acc_beta_bdf[eN_k],
		     mom_w_acc,
		     dmom_w_acc_w,
		     mom_w_acc_t,
		     dmom_w_acc_w_t);
	      //
	      //calculate subgrid error (strong residual and adjoint)
	      //
	      //calculate strong residual
	      pdeResidual_p = ck.Mass_strong(mass_acc_t) +
		ck.Advection_strong(dmass_adv_u,grad_u) +
		ck.Advection_strong(dmass_adv_v,grad_v) +
		ck.Advection_strong(dmass_adv_w,grad_w);
	  
              dmom_adv_sge[0] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+0] - MOVING_DOMAIN*xt);
              dmom_adv_sge[1] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+1] - MOVING_DOMAIN*yt);
              dmom_adv_sge[2] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+2] - MOVING_DOMAIN*zt);

	      pdeResidual_u = ck.Mass_strong(mom_u_acc_t) +
		ck.Advection_strong(dmom_adv_sge,grad_u) +
		ck.Hamiltonian_strong(dmom_u_ham_grad_p,grad_p) +
		ck.Reaction_strong(mom_u_source);
	  
	      pdeResidual_v = ck.Mass_strong(mom_v_acc_t) +
		ck.Advection_strong(dmom_adv_sge,grad_v) +
		ck.Hamiltonian_strong(dmom_v_ham_grad_p,grad_p) + 
		ck.Reaction_strong(mom_v_source);
	  
	      pdeResidual_w = ck.Mass_strong(mom_w_acc_t) + 
		ck.Advection_strong(dmom_adv_sge,grad_w) +
		ck.Hamiltonian_strong(dmom_w_ham_grad_p,grad_p) +
		ck.Reaction_strong(mom_w_source);
	
	      //calculate tau and tau*Res
	      //cek check whether we need to change tau based on new terms
	      calculateSubgridError_tau(hFactor,
					elementDiameter[eN],
					dmom_u_acc_u_t,
					dmom_u_acc_u,
					dmom_adv_sge,
					mom_u_diff_ten[1],
					dmom_u_ham_grad_p[0],
					tau_v0,
					tau_p0,
					q_cfl[eN_k]);


	      calculateSubgridError_tau(Ct_sge,Cd_sge,
			                G,G_dd_G,tr_G,
					dmom_u_acc_u_t,
					dmom_adv_sge,
					mom_u_diff_ten[1],
					dmom_u_ham_grad_p[0],
					tau_v1,
					tau_p1,
					q_cfl[eN_k]);	


	      tau_v = useMetrics*tau_v1+(1.0-useMetrics)*tau_v0;
	      tau_p = useMetrics*tau_p1+(1.0-useMetrics)*tau_p0;


	      //calculateSubgridError_tau(Ct_sge,
	      // 				      Cd_sge,
	      // 				      G,
	      // 				      G_dd_G,
	      // 				      tr_G,
	      // 				      dmom_u_acc_u,//rho
	      // 				      dmom_u_acc_u_t/dmom_u_acc_u,//Dt
	      // 				      &q_velocity_sge[eN_k_nSpace],//v
	      // 				      mom_u_diff_ten[1],//mu
	      // 				      tau_p,
	      // 				      tau_v,
	      // 				      q_cfl[eN_k]);
	      calculateSubgridError_tauRes(tau_p,
					   tau_v,
					   pdeResidual_p,
					   pdeResidual_u,
					   pdeResidual_v,
					   pdeResidual_w,
					   subgridError_p,
					   subgridError_u,
					   subgridError_v,
					   subgridError_w);
	      // velocity used in adjoint (VMS or RBLES, with or without lagging the grid scale velocity)
	      dmom_adv_star[0] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+0] - MOVING_DOMAIN*xt + useRBLES*subgridError_u);
	      dmom_adv_star[1] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+1] - MOVING_DOMAIN*yt + useRBLES*subgridError_v);
              dmom_adv_star[2] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+2] - MOVING_DOMAIN*zt + useRBLES*subgridError_w);
        
	      // adjoint times the test functions 
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int i_nSpace = i*nSpace;
		  Lstar_u_p[i]=ck.Advection_adjoint(dmass_adv_u,&p_grad_test_dV[i_nSpace]);
		  Lstar_v_p[i]=ck.Advection_adjoint(dmass_adv_v,&p_grad_test_dV[i_nSpace]);
		  Lstar_w_p[i]=ck.Advection_adjoint(dmass_adv_w,&p_grad_test_dV[i_nSpace]);
                  //use the same advection adjoint for all three since we're approximating the linearized adjoint
		  Lstar_u_u[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
		  Lstar_v_v[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
		  Lstar_w_w[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
		  Lstar_p_u[i]=ck.Hamiltonian_adjoint(dmom_u_ham_grad_p,&vel_grad_test_dV[i_nSpace]) + ck.Reaction_adjoint(dmom_u_source_p,vel_test_dV[i]);
		  Lstar_p_v[i]=ck.Hamiltonian_adjoint(dmom_v_ham_grad_p,&vel_grad_test_dV[i_nSpace]) + ck.Reaction_adjoint(dmom_v_source_p,vel_test_dV[i]);
		  Lstar_p_w[i]=ck.Hamiltonian_adjoint(dmom_w_ham_grad_p,&vel_grad_test_dV[i_nSpace]) + ck.Reaction_adjoint(dmom_w_source_p,vel_test_dV[i]);
		}

	      norm_Rv = sqrt(pdeResidual_u*pdeResidual_u + pdeResidual_v*pdeResidual_v + pdeResidual_w*pdeResidual_w);
	      q_numDiff_u[eN_k] = C_dc*norm_Rv*(useMetrics/sqrt(G_dd_G)  + 
	                                        (1.0-useMetrics)*hFactor*hFactor*elementDiameter[eN]*elementDiameter[eN]);
	      q_numDiff_v[eN_k] = q_numDiff_u[eN_k];
	      q_numDiff_w[eN_k] = q_numDiff_u[eN_k];
	      // //cek debug
	      // q_numDiff_u[eN_k] = 0.0;
	      // q_numDiff_v[eN_k] = 0.0;
	      // q_numDiff_w[eN_k] = 0.0;
	      // q_numDiff_u_last[eN_k] = 0.0;	  
	      // q_numDiff_v_last[eN_k] = 0.0;	  
	      // q_numDiff_w_last[eN_k] = 0.0;	  
	      // 
	      //update element residual 
	      // 
	      for(int i=0;i<nDOF_test_element;i++) 
		{ 
		  register int i_nSpace=i*nSpace;

		  elementResidual_p[i] += ck.Mass_weak(mass_acc_t,p_test_dV[i]) +
		    ck.Advection_weak(mass_adv,&p_grad_test_dV[i_nSpace]) +
		    ck.SubgridError(subgridError_u,Lstar_u_p[i]) + 
		    ck.SubgridError(subgridError_v,Lstar_v_p[i]) + 
		    ck.SubgridError(subgridError_w,Lstar_w_p[i]);

		  elementResidual_u[i] += ck.Mass_weak(mom_u_acc_t,vel_test_dV[i]) + 
		    ck.Advection_weak(mom_u_adv,&vel_grad_test_dV[i_nSpace]) +
		    ck.Diffusion_weak(sdInfo_u_u_rowptr,sdInfo_u_u_colind,mom_u_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Diffusion_weak(sdInfo_u_v_rowptr,sdInfo_u_v_colind,mom_uv_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Diffusion_weak(sdInfo_u_w_rowptr,sdInfo_u_w_colind,mom_uw_diff_ten,grad_w,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Reaction_weak(mom_u_source,vel_test_dV[i]) + 
		    ck.Hamiltonian_weak(mom_u_ham,vel_test_dV[i]) + 
		    ck.SubgridError(subgridError_p,Lstar_p_u[i]) + 
		    ck.SubgridError(subgridError_u,Lstar_u_u[i]) + 
		    ck.NumericalDiffusion(q_numDiff_u_last[eN_k],grad_u,&vel_grad_test_dV[i_nSpace]); 
		 
		  elementResidual_v[i] += ck.Mass_weak(mom_v_acc_t,vel_test_dV[i]) + 
		    ck.Advection_weak(mom_v_adv,&vel_grad_test_dV[i_nSpace]) +
		    ck.Diffusion_weak(sdInfo_v_v_rowptr,sdInfo_v_v_colind,mom_v_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Diffusion_weak(sdInfo_v_u_rowptr,sdInfo_v_u_colind,mom_vu_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Diffusion_weak(sdInfo_v_w_rowptr,sdInfo_v_w_colind,mom_vw_diff_ten,grad_w,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Reaction_weak(mom_v_source,vel_test_dV[i]) + 
		    ck.Hamiltonian_weak(mom_v_ham,vel_test_dV[i]) + 
		    ck.SubgridError(subgridError_p,Lstar_p_v[i]) + 
		    ck.SubgridError(subgridError_v,Lstar_v_v[i]) + 
		    ck.NumericalDiffusion(q_numDiff_v_last[eN_k],grad_v,&vel_grad_test_dV[i_nSpace]); 

		  elementResidual_w[i] +=  ck.Mass_weak(mom_w_acc_t,vel_test_dV[i]) +
		    ck.Advection_weak(mom_w_adv,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Diffusion_weak(sdInfo_w_w_rowptr,sdInfo_w_w_colind,mom_w_diff_ten,grad_w,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Diffusion_weak(sdInfo_w_u_rowptr,sdInfo_w_u_colind,mom_wu_diff_ten,grad_u,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Diffusion_weak(sdInfo_w_v_rowptr,sdInfo_w_v_colind,mom_wv_diff_ten,grad_v,&vel_grad_test_dV[i_nSpace]) + 
		    ck.Reaction_weak(mom_w_source,vel_test_dV[i]) + 
		    ck.Hamiltonian_weak(mom_w_ham,vel_test_dV[i]) + 
		    ck.SubgridError(subgridError_p,Lstar_p_w[i]) + 
		    ck.SubgridError(subgridError_w,Lstar_w_w[i]) + 
		    ck.NumericalDiffusion(q_numDiff_w_last[eN_k],grad_w,&vel_grad_test_dV[i_nSpace]); 
		}//i
	    }
	  //
	  //load element into global residual and save element residual
	  //
	  for(int i=0;i<nDOF_test_element;i++) 
	    { 
	      register int eN_i=eN*nDOF_test_element+i;

	      elementResidual_p_save[eN_i] +=  elementResidual_p[i];
	  
	      globalResidual[offset_p+stride_p*p_l2g[eN_i]]+=elementResidual_p[i];
	      globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
	      globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
	      globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
	    }//i
	  // //
	  // //debug
	  // //
	  // for(int i=0;i<nDOF_test_element;i++) 
	  //   { 
	  // 	  std::cout<<"eN "<<eN<<" i "<<i<<std::endl;
	  // 	  std::cout<<"r_p"<<elementResidual_p[i]<<std::endl;
	  // 	  std::cout<<"r_u"<<elementResidual_u[i]<<std::endl;
	  // 	  std::cout<<"r_v"<<elementResidual_v[i]<<std::endl;
	  // 	  std::cout<<"r_w"<<elementResidual_w[i]<<std::endl;
	  // 	}
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
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element;
	  register double elementResidual_p[nDOF_test_element],
	    elementResidual_u[nDOF_test_element],
	    elementResidual_v[nDOF_test_element],
	    elementResidual_w[nDOF_test_element],
	    eps_rho,eps_mu;
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      elementResidual_p[i]=0.0;
	      elementResidual_u[i]=0.0;
	      elementResidual_v[i]=0.0;
	      elementResidual_w[i]=0.0;
	    }
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;
	      register double p_ext=0.0,C_ext=0.0,
		u_ext=0.0,
		v_ext=0.0,
		w_ext=0.0,
		grad_p_ext[nSpace],
		grad_C_ext[nSpace],
		grad_u_ext[nSpace],
		grad_v_ext[nSpace],
		grad_w_ext[nSpace],
		mass_acc_ext=0.0,
		dmass_acc_ext=0.0,
		mom_u_acc_ext=0.0,
		dmom_u_acc_u_ext=0.0,
		mom_v_acc_ext=0.0,
		dmom_v_acc_v_ext=0.0,
		mom_w_acc_ext=0.0,
		dmom_w_acc_w_ext=0.0,
		mass_adv_ext[nSpace],
		dmass_adv_u_ext[nSpace],
		dmass_adv_v_ext[nSpace],
		dmass_adv_w_ext[nSpace],
		mom_u_adv_ext[nSpace],
		dmom_u_adv_u_ext[nSpace],
		dmom_u_adv_v_ext[nSpace],
		dmom_u_adv_w_ext[nSpace],
		mom_v_adv_ext[nSpace],
		dmom_v_adv_u_ext[nSpace],
		dmom_v_adv_v_ext[nSpace],
		dmom_v_adv_w_ext[nSpace],
		mom_w_adv_ext[nSpace],
		dmom_w_adv_u_ext[nSpace],
		dmom_w_adv_v_ext[nSpace],
		dmom_w_adv_w_ext[nSpace],
		mom_u_diff_ten_ext[nSpace],
		mom_v_diff_ten_ext[nSpace],
		mom_w_diff_ten_ext[nSpace],
		mom_uv_diff_ten_ext[1],
		mom_uw_diff_ten_ext[1],
		mom_vu_diff_ten_ext[1],
		mom_vw_diff_ten_ext[1],
		mom_wu_diff_ten_ext[1],
		mom_wv_diff_ten_ext[1],
		mom_u_source_ext=0.0,
		dmom_u_source_p_ext=0.0,
		mom_v_source_ext=0.0,
		dmom_v_source_p_ext=0.0,
		mom_w_source_ext=0.0,
		dmom_w_source_p_ext=0.0,
		mom_u_ham_ext=0.0,
		dmom_u_ham_grad_p_ext[nSpace],
		mom_v_ham_ext=0.0,
		dmom_v_ham_grad_p_ext[nSpace],
		mom_w_ham_ext=0.0,
		dmom_w_ham_grad_p_ext[nSpace],
		dmom_u_adv_p_ext[nSpace],
		dmom_v_adv_p_ext[nSpace],
		dmom_w_adv_p_ext[nSpace],
		flux_mass_ext=0.0,
		flux_mom_u_adv_ext=0.0,
		flux_mom_v_adv_ext=0.0,
		flux_mom_w_adv_ext=0.0,
		flux_mom_u_diff_ext=0.0,
		flux_mom_v_diff_ext=0.0,
		flux_mom_w_diff_ext=0.0,
		bc_p_ext=0.0,
		bc_C_ext=0.0,
		bc_u_ext=0.0,
		bc_v_ext=0.0,
		bc_w_ext=0.0,
		bc_mass_acc_ext=0.0,
		bc_dmass_acc_ext=0.0,
		bc_mom_u_acc_ext=0.0,
		bc_dmom_u_acc_u_ext=0.0,
		bc_mom_v_acc_ext=0.0,
		bc_dmom_v_acc_v_ext=0.0,
		bc_mom_w_acc_ext=0.0,
		bc_dmom_w_acc_w_ext=0.0,
		bc_mass_adv_ext[nSpace],
		bc_dmass_adv_u_ext[nSpace],
		bc_dmass_adv_v_ext[nSpace],
		bc_dmass_adv_w_ext[nSpace],
		bc_mom_u_adv_ext[nSpace],
		bc_dmom_u_adv_u_ext[nSpace],
		bc_dmom_u_adv_v_ext[nSpace],
		bc_dmom_u_adv_w_ext[nSpace],
		bc_mom_v_adv_ext[nSpace],
		bc_dmom_v_adv_u_ext[nSpace],
		bc_dmom_v_adv_v_ext[nSpace],
		bc_dmom_v_adv_w_ext[nSpace],
		bc_mom_w_adv_ext[nSpace],
		bc_dmom_w_adv_u_ext[nSpace],
		bc_dmom_w_adv_v_ext[nSpace],
		bc_dmom_w_adv_w_ext[nSpace],
		bc_mom_u_diff_ten_ext[nSpace],
		bc_mom_v_diff_ten_ext[nSpace],
		bc_mom_w_diff_ten_ext[nSpace],
		bc_mom_uv_diff_ten_ext[1],
		bc_mom_uw_diff_ten_ext[1],
		bc_mom_vu_diff_ten_ext[1],
		bc_mom_vw_diff_ten_ext[1],
		bc_mom_wu_diff_ten_ext[1],
		bc_mom_wv_diff_ten_ext[1],
		bc_mom_u_source_ext=0.0,
		bc_mom_v_source_ext=0.0,
		bc_mom_w_source_ext=0.0,
		bc_mom_u_ham_ext=0.0,
		bc_dmom_u_ham_grad_p_ext[nSpace],
		bc_mom_v_ham_ext=0.0,
		bc_dmom_v_ham_grad_p_ext[nSpace],
		bc_mom_w_ham_ext=0.0,
		bc_dmom_w_ham_grad_p_ext[nSpace],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		dS,p_test_dS[nDOF_test_element],vel_test_dS[nDOF_test_element],
		p_grad_trial_trace[nDOF_trial_element*nSpace],vel_grad_trial_trace[nDOF_trial_element*nSpace],
		normal[3],x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty;
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
	      //xt_ext=0.0;yt_ext=0.0;zt_ext=0.0;
	      //std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<zt_ext<<std::endl;
	      //std::cout<<"integralScaling - metricTensorDetSrt ==============================="<<integralScaling-metricTensorDetSqrt<<std::endl;
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      //get the metric tensor
	      //cek todo use symmetry
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      ck.calculateGScale(G,&ebqe_normal_phi_ext[ebNE_kb_nSpace],h_phi);
	      
	      eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      
	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&p_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
	      ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],p_ext);
	      ck.valFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],C_ext);
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],w_ext);
	      ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
	      ck.gradFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_C_ext);
	      ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_u_ext);
	      ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_v_ext);
	      ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_w_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  p_test_dS[j] = p_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		}
	      // //
	      // //debugging section for finite element calculations on exterior
	      // //
	      // std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl;
	      // for (int j=0;j<nDOF_trial_element;j++)
	      //   {
	      //     std::cout<<"p_trial_trace["<<j<<"] "<<p_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	      //     std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	      //     std::cout<<"p_test_dS["<<j<<"] "<<p_test_dS[j]<<std::endl;
	      //     std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl;
	      //     for (int I=0;I<nSpace;I++)
	      // 	{
	      // 	  std::cout<<"p_grad_trial_trace["<<j<<","<<I<<"] "<<p_grad_trial_trace[j*nSpace+I]<<std::endl;
	      // 	  std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl;
	      // 	}
	      //   }
	      // std::cout<<"p_ext "<<p_ext<<std::endl;
	      // std::cout<<"u_ext "<<u_ext<<std::endl;
	      // std::cout<<"v_ext "<<v_ext<<std::endl;
	      // std::cout<<"w_ext "<<w_ext<<std::endl;
	      // for(int I=0;I<nSpace;I++)
	      //   {
	      //     std::cout<<"grad_p_ext["<<I<<"] "<<grad_p_ext[I]<<std::endl;
	      //     std::cout<<"grad_u_ext["<<I<<"] "<<grad_u_ext[I]<<std::endl;
	      //     std::cout<<"grad_v_ext["<<I<<"] "<<grad_v_ext[I]<<std::endl;
	      //     std::cout<<"grad_w_ext["<<I<<"] "<<grad_w_ext[I]<<std::endl;
	      //   }
	      //
	      //load the boundary values
	      //
	      bc_p_ext = isDOFBoundary_p[ebNE_kb]*ebqe_bc_p_ext[ebNE_kb]+(1-isDOFBoundary_p[ebNE_kb])*p_ext;
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	      bc_v_ext = isDOFBoundary_v[ebNE_kb]*ebqe_bc_v_ext[ebNE_kb]+(1-isDOFBoundary_v[ebNE_kb])*v_ext;
	      bc_w_ext = isDOFBoundary_w[ebNE_kb]*ebqe_bc_w_ext[ebNE_kb]+(1-isDOFBoundary_w[ebNE_kb])*w_ext;
	      // 
	      //calculate the pde coefficients using the solution and the boundary values for the solution 
	      // 
	      //cek debug
	      //eps_rho=0.1;
	      //eps_mu=0.1;
	      double ebqe_sedimentForce[nSpace];
	      evaluateCoefficients(eps_rho,
				   eps_mu,
				   sigma,
				   rho_0,
				   nu_0,
				   rho_1,
				   nu_1,
				   g,
				   ebqe_sedimentForce,//[ebNE_kb_nSpace],
				   ebqe_phi_ext[ebNE_kb],
				   &ebqe_normal_phi_ext[ebNE_kb_nSpace],
				   ebqe_kappa_phi_ext[ebNE_kb],
				   C_ext,
				   grad_C_ext,
				   p_ext,
				   grad_p_ext,
				   u_ext,
				   v_ext,
				   w_ext,
				   mass_acc_ext,
				   mom_u_acc_ext,
				   dmom_u_acc_u_ext,
				   mom_v_acc_ext,
				   dmom_v_acc_v_ext,
				   mom_w_acc_ext,
				   dmom_w_acc_w_ext,
				   mass_adv_ext,
				   dmass_adv_u_ext,
				   dmass_adv_v_ext,
				   dmass_adv_w_ext,
				   mom_u_adv_ext,
				   dmom_u_adv_u_ext,
				   dmom_u_adv_v_ext,
				   dmom_u_adv_w_ext,
				   mom_v_adv_ext,
				   dmom_v_adv_u_ext,
				   dmom_v_adv_v_ext,
				   dmom_v_adv_w_ext,
				   mom_w_adv_ext,
				   dmom_w_adv_u_ext,
				   dmom_w_adv_v_ext,
				   dmom_w_adv_w_ext,
				   mom_u_diff_ten_ext,
				   mom_v_diff_ten_ext,
				   mom_w_diff_ten_ext,
				   mom_uv_diff_ten_ext,
				   mom_uw_diff_ten_ext,
				   mom_vu_diff_ten_ext,
				   mom_vw_diff_ten_ext,
				   mom_wu_diff_ten_ext,
				   mom_wv_diff_ten_ext,
				   mom_u_source_ext,
				   dmom_u_source_p_ext,
				   mom_v_source_ext,
				   dmom_v_source_p_ext,
				   mom_w_source_ext,
				   dmom_w_source_p_ext,
				   mom_u_ham_ext,
				   dmom_u_ham_grad_p_ext,
				   mom_v_ham_ext,
				   dmom_v_ham_grad_p_ext,
				   mom_w_ham_ext,
				   dmom_w_ham_grad_p_ext);          
	      evaluateCoefficients(eps_rho,
				   eps_mu,
				   sigma,
				   rho_0,
				   nu_0,
				   rho_1,
				   nu_1,
				   g,
				   ebqe_sedimentForce,//[ebNE_kb_nSpace],
				   ebqe_phi_ext[ebNE_kb],
				   &ebqe_normal_phi_ext[ebNE_kb_nSpace],
				   ebqe_kappa_phi_ext[ebNE_kb],
				   bc_C_ext,
				   grad_C_ext,//cek shouldn't be used
				   bc_p_ext,
				   grad_p_ext,//cek should't be used
				   bc_u_ext,
				   bc_v_ext,
				   bc_w_ext,
				   bc_mass_acc_ext,
				   bc_mom_u_acc_ext,
				   bc_dmom_u_acc_u_ext,
				   bc_mom_v_acc_ext,
				   bc_dmom_v_acc_v_ext,
				   bc_mom_w_acc_ext,
				   bc_dmom_w_acc_w_ext,
				   bc_mass_adv_ext,
				   bc_dmass_adv_u_ext,
				   bc_dmass_adv_v_ext,
				   bc_dmass_adv_w_ext,
				   bc_mom_u_adv_ext,
				   bc_dmom_u_adv_u_ext,
				   bc_dmom_u_adv_v_ext,
				   bc_dmom_u_adv_w_ext,
				   bc_mom_v_adv_ext,
				   bc_dmom_v_adv_u_ext,
				   bc_dmom_v_adv_v_ext,
				   bc_dmom_v_adv_w_ext,
				   bc_mom_w_adv_ext,
				   bc_dmom_w_adv_u_ext,
				   bc_dmom_w_adv_v_ext,
				   bc_dmom_w_adv_w_ext,
				   bc_mom_u_diff_ten_ext,
				   bc_mom_v_diff_ten_ext,
				   bc_mom_w_diff_ten_ext,
				   bc_mom_uv_diff_ten_ext,
				   bc_mom_uw_diff_ten_ext,
				   bc_mom_vu_diff_ten_ext,
				   bc_mom_vw_diff_ten_ext,
				   bc_mom_wu_diff_ten_ext,
				   bc_mom_wv_diff_ten_ext,
				   bc_mom_u_source_ext,
				   dmom_u_source_p_ext,//cek not used
				   bc_mom_v_source_ext,
				   dmom_v_source_p_ext,//cek not used
				   bc_mom_w_source_ext,
				   dmom_w_source_p_ext,//cek not used
				   bc_mom_u_ham_ext,
				   bc_dmom_u_ham_grad_p_ext,
				   bc_mom_v_ham_ext,
				   bc_dmom_v_ham_grad_p_ext,
				   bc_mom_w_ham_ext,
				   bc_dmom_w_ham_grad_p_ext);          
	      //
	      //moving domain
	      //
	      mass_adv_ext[0] -= MOVING_DOMAIN*mass_acc_ext*xt_ext;
	      mass_adv_ext[1] -= MOVING_DOMAIN*mass_acc_ext*yt_ext;
	      mass_adv_ext[2] -= MOVING_DOMAIN*mass_acc_ext*zt_ext;

	      mom_u_adv_ext[0] -= MOVING_DOMAIN*mom_u_acc_ext*xt_ext;
	      mom_u_adv_ext[1] -= MOVING_DOMAIN*mom_u_acc_ext*yt_ext;
	      mom_u_adv_ext[2] -= MOVING_DOMAIN*mom_u_acc_ext*zt_ext;
	      dmom_u_adv_u_ext[0] -= MOVING_DOMAIN*dmom_u_acc_u_ext*xt_ext;
	      dmom_u_adv_u_ext[1] -= MOVING_DOMAIN*dmom_u_acc_u_ext*yt_ext;
	      dmom_u_adv_u_ext[2] -= MOVING_DOMAIN*dmom_u_acc_u_ext*zt_ext;

	      mom_v_adv_ext[0] -= MOVING_DOMAIN*mom_v_acc_ext*xt_ext;
	      mom_v_adv_ext[1] -= MOVING_DOMAIN*mom_v_acc_ext*yt_ext;
	      mom_v_adv_ext[2] -= MOVING_DOMAIN*mom_v_acc_ext*zt_ext;
	      dmom_v_adv_v_ext[0] -= MOVING_DOMAIN*dmom_v_acc_v_ext*xt_ext;
	      dmom_v_adv_v_ext[1] -= MOVING_DOMAIN*dmom_v_acc_v_ext*yt_ext;
	      dmom_v_adv_v_ext[2] -= MOVING_DOMAIN*dmom_v_acc_v_ext*zt_ext;

	      mom_w_adv_ext[0] -= MOVING_DOMAIN*mom_w_acc_ext*xt_ext;
	      mom_w_adv_ext[1] -= MOVING_DOMAIN*mom_w_acc_ext*yt_ext;
	      mom_w_adv_ext[2] -= MOVING_DOMAIN*mom_w_acc_ext*zt_ext;
	      dmom_w_adv_w_ext[0] -= MOVING_DOMAIN*dmom_w_acc_w_ext*xt_ext;
	      dmom_w_adv_w_ext[1] -= MOVING_DOMAIN*dmom_w_acc_w_ext*yt_ext;
	      dmom_w_adv_w_ext[2] -= MOVING_DOMAIN*dmom_w_acc_w_ext*zt_ext;

	      //bc's
	      bc_mom_u_adv_ext[0] -= MOVING_DOMAIN*bc_mom_u_acc_ext*xt_ext;
	      bc_mom_u_adv_ext[1] -= MOVING_DOMAIN*bc_mom_u_acc_ext*yt_ext;
	      bc_mom_u_adv_ext[2] -= MOVING_DOMAIN*bc_mom_u_acc_ext*zt_ext;

	      bc_mom_v_adv_ext[0] -= MOVING_DOMAIN*bc_mom_v_acc_ext*xt_ext;
	      bc_mom_v_adv_ext[1] -= MOVING_DOMAIN*bc_mom_v_acc_ext*yt_ext;
	      bc_mom_v_adv_ext[2] -= MOVING_DOMAIN*bc_mom_v_acc_ext*zt_ext;

	      bc_mom_w_adv_ext[0] -= MOVING_DOMAIN*bc_mom_w_acc_ext*xt_ext;
	      bc_mom_w_adv_ext[1] -= MOVING_DOMAIN*bc_mom_w_acc_ext*yt_ext;
	      bc_mom_w_adv_ext[2] -= MOVING_DOMAIN*bc_mom_w_acc_ext*zt_ext;
	      // 
	      //calculate the numerical fluxes 
	      // 
	      //cek debug
	      //ebqe_penalty_ext[ebNE_kb] = 10.0;
	      //
	      ck.calculateGScale(G,normal,h_penalty);
	      h_penalty = 10.0/h_penalty;
	      //cek debug, do it the old way
	      h_penalty = 100.0/elementDiameter[eN];
	      exteriorNumericalAdvectiveFlux(isDOFBoundary_p[ebNE_kb],
					     isDOFBoundary_u[ebNE_kb],
					     isDOFBoundary_v[ebNE_kb],
					     isDOFBoundary_w[ebNE_kb],
					     isAdvectiveFluxBoundary_p[ebNE_kb],
					     isAdvectiveFluxBoundary_u[ebNE_kb],
					     isAdvectiveFluxBoundary_v[ebNE_kb],
					     isAdvectiveFluxBoundary_w[ebNE_kb],
					     dmom_u_ham_grad_p_ext[0],//=1/rho,
					     normal,
					     bc_p_ext,
					     bc_mass_adv_ext,
					     bc_mom_u_adv_ext,
					     bc_mom_v_adv_ext,
					     bc_mom_w_adv_ext,
					     ebqe_bc_flux_mass_ext[ebNE_kb],
					     ebqe_bc_flux_mom_u_adv_ext[ebNE_kb],
					     ebqe_bc_flux_mom_v_adv_ext[ebNE_kb],
					     ebqe_bc_flux_mom_w_adv_ext[ebNE_kb],
					     p_ext,
					     mass_adv_ext,
					     mom_u_adv_ext,
					     mom_v_adv_ext,
					     mom_w_adv_ext,
					     dmass_adv_u_ext,
					     dmass_adv_v_ext,
					     dmass_adv_w_ext,
					     dmom_u_adv_p_ext,
					     dmom_u_adv_u_ext,
					     dmom_u_adv_v_ext,
					     dmom_u_adv_w_ext,
					     dmom_v_adv_p_ext,
					     dmom_v_adv_u_ext,
					     dmom_v_adv_v_ext,
					     dmom_v_adv_w_ext,
					     dmom_w_adv_p_ext,
					     dmom_w_adv_u_ext,
					     dmom_w_adv_v_ext,
					     dmom_w_adv_w_ext,
					     flux_mass_ext,
					     flux_mom_u_adv_ext,
					     flux_mom_v_adv_ext,
					     flux_mom_w_adv_ext,
					     &ebqe_velocity[ebNE_kb_nSpace]);
	      //cek todo need to switch to full stress and add adjoint consistency
	      exteriorNumericalDiffusiveFlux(eps_rho,
					     ebqe_phi_ext[ebNE_kb],
					     sdInfo_u_u_rowptr,
					     sdInfo_u_u_colind,
					     isDOFBoundary_u[ebNE_kb],
					     isDiffusiveFluxBoundary_u[ebNE_kb],
					     normal,
					     bc_mom_u_diff_ten_ext,
					     bc_u_ext,
					     ebqe_bc_flux_u_diff_ext[ebNE_kb],
					     mom_u_diff_ten_ext,
					     grad_u_ext,
					     u_ext,
					     h_penalty,//ebqe_penalty_ext[ebNE_kb],
					     flux_mom_u_diff_ext);
	      exteriorNumericalDiffusiveFlux(eps_rho,
					     ebqe_phi_ext[ebNE_kb],
					     sdInfo_v_v_rowptr,
					     sdInfo_v_v_colind,
					     isDOFBoundary_v[ebNE_kb],
					     isDiffusiveFluxBoundary_v[ebNE_kb],
					     normal,
					     bc_mom_v_diff_ten_ext,
					     bc_v_ext,
					     ebqe_bc_flux_v_diff_ext[ebNE_kb],
					     mom_v_diff_ten_ext,
					     grad_v_ext,
					     v_ext,
					     h_penalty,//ebqe_penalty_ext[ebNE_kb],
					     flux_mom_v_diff_ext);
	      exteriorNumericalDiffusiveFlux(eps_rho,
					     ebqe_phi_ext[ebNE_kb],
					     sdInfo_w_w_rowptr,
					     sdInfo_w_w_colind,
					     isDOFBoundary_w[ebNE_kb],
					     isDiffusiveFluxBoundary_w[ebNE_kb],
					     normal,
					     bc_mom_w_diff_ten_ext,
					     bc_w_ext,
					     ebqe_bc_flux_w_diff_ext[ebNE_kb],
					     mom_w_diff_ten_ext,
					     grad_w_ext,
					     w_ext,
					     h_penalty,//ebqe_penalty_ext[ebNE_kb],
					     flux_mom_w_diff_ext);
	      flux[ebN*nQuadraturePoints_elementBoundary+kb] = flux_mass_ext;
	      //flux[ebN*nQuadraturePoints_elementBoundary+kb] = 0.0;//cek debug
	      //
	      //update residuals
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  elementResidual_p[i] += ck.ExteriorElementBoundaryFlux(flux_mass_ext,p_test_dS[i]);
		  globalConservationError += ck.ExteriorElementBoundaryFlux(flux_mass_ext,p_test_dS[i]);
		  
		  elementResidual_u[i] += ck.ExteriorElementBoundaryFlux(flux_mom_u_adv_ext,vel_test_dS[i])+
		    ck.ExteriorElementBoundaryFlux(flux_mom_u_diff_ext,vel_test_dS[i]); 

		  elementResidual_v[i] += ck.ExteriorElementBoundaryFlux(flux_mom_v_adv_ext,vel_test_dS[i]) +
		    ck.ExteriorElementBoundaryFlux(flux_mom_v_diff_ext,vel_test_dS[i]); 
	       
		  elementResidual_w[i] += ck.ExteriorElementBoundaryFlux(flux_mom_w_adv_ext,vel_test_dS[i]) +
		    ck.ExteriorElementBoundaryFlux(flux_mom_w_diff_ext,vel_test_dS[i]); 
		}//i
	    }//kb
	  //
	  //update the element and global residual storage
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      int eN_i = eN*nDOF_test_element+i;
	      
	      elementResidual_p_save[eN_i] +=  elementResidual_p[i];
		  
	      globalResidual[offset_p+stride_p*p_l2g[eN_i]]+=elementResidual_p[i];
	      globalResidual[offset_u+stride_u*vel_l2g[eN_i]]+=elementResidual_u[i];
	      globalResidual[offset_v+stride_v*vel_l2g[eN_i]]+=elementResidual_v[i];
	      globalResidual[offset_w+stride_w*vel_l2g[eN_i]]+=elementResidual_w[i];
	    }//i
	  // //
	  // //debug
	  // //
	  // for(int i=0;i<nDOF_test_element;i++) 
	  //   { 
	  // 	  std::cout<<"ebNE "<<ebNE<<" i "<<i<<std::endl;
	  // 	  std::cout<<"r_p"<<elementResidual_p[i]<<std::endl;
	  // 	  std::cout<<"r_u"<<elementResidual_u[i]<<std::endl;
	  // 	  std::cout<<"r_v"<<elementResidual_v[i]<<std::endl;
	  // 	  std::cout<<"r_w"<<elementResidual_w[i]<<std::endl;
	  // 	}

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
			   double* p_trial_ref,
			   double* p_grad_trial_ref,
			   double* p_test_ref,
			   double* p_grad_test_ref,
			   double* vel_trial_ref,
			   double* vel_grad_trial_ref,
			   double* vel_test_ref,
			   double* vel_grad_test_ref,
			   //element boundary
			   double* mesh_trial_trace_ref,
			   double* mesh_grad_trial_trace_ref,
			   double* dS_ref,
			   double* p_trial_trace_ref,
			   double* p_grad_trial_trace_ref,
			   double* p_test_trace_ref,
			   double* p_grad_test_trace_ref,
			   double* vel_trial_trace_ref,
			   double* vel_grad_trial_trace_ref,
			   double* vel_test_trace_ref,
			   double* vel_grad_test_trace_ref,					 
			   double* normal_ref,
			   double* boundaryJac_ref,
			   //physics
			   double* elementDiameter,
			   double hFactor,
			   int nElements_global,
			   double useRBLES,
			   double useMetrics, 
			   double alphaBDF,
			   double epsFact_rho,
			   double epsFact_mu,
			   double sigma,
			   double rho_0,
			   double nu_0,
			   double rho_1,
			   double nu_1,
			   double Ct_sge,
			   double Cd_sge,
			   double C_dg,
			   int* p_l2g, 
			   int* vel_l2g,
			   double* C_dof, double* p_dof, double* u_dof, double* v_dof, double* w_dof,
			   double* g,
			   double* q_sedimentForce,
			   //double* ebqe_sedimentForce,
			   double* phi,
			   double* normal_phi,
			   double* kappa_phi,
			   double* q_mass_acc_beta_bdf, double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
			   double* q_velocity_sge,
			   double* q_cfl,
			   double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
			   int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
			   int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
			   int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
			   int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
			   int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
			   int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
			   int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
			   int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
			   int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
			   int* csrRowIndeces_p_p,int* csrColumnOffsets_p_p,
			   int* csrRowIndeces_p_u,int* csrColumnOffsets_p_u,
			   int* csrRowIndeces_p_v,int* csrColumnOffsets_p_v,
			   int* csrRowIndeces_p_w,int* csrColumnOffsets_p_w,
			   int* csrRowIndeces_u_p,int* csrColumnOffsets_u_p,
			   int* csrRowIndeces_u_u,int* csrColumnOffsets_u_u,
			   int* csrRowIndeces_u_v,int* csrColumnOffsets_u_v,
			   int* csrRowIndeces_u_w,int* csrColumnOffsets_u_w,
			   int* csrRowIndeces_v_p,int* csrColumnOffsets_v_p,
			   int* csrRowIndeces_v_u,int* csrColumnOffsets_v_u,
			   int* csrRowIndeces_v_v,int* csrColumnOffsets_v_v,
			   int* csrRowIndeces_v_w,int* csrColumnOffsets_v_w,
			   int* csrRowIndeces_w_p,int* csrColumnOffsets_w_p,
			   int* csrRowIndeces_w_u,int* csrColumnOffsets_w_u,
			   int* csrRowIndeces_w_v,int* csrColumnOffsets_w_v,
			   int* csrRowIndeces_w_w,int* csrColumnOffsets_w_w,
			   double* globalJacobian,
			   int nExteriorElementBoundaries_global,
			   int* exteriorElementBoundariesArray,
			   int* elementBoundaryElementsArray,
			   int* elementBoundaryLocalElementBoundariesArray,
			   double* ebqe_phi_ext,
			   double* ebqe_normal_phi_ext,
			   double* ebqe_kappa_phi_ext,
			   int* isDOFBoundary_p,
			   int* isDOFBoundary_u,
			   int* isDOFBoundary_v,
			   int* isDOFBoundary_w,
			   int* isAdvectiveFluxBoundary_p,
			   int* isAdvectiveFluxBoundary_u,
			   int* isAdvectiveFluxBoundary_v,
			   int* isAdvectiveFluxBoundary_w,
			   int* isDiffusiveFluxBoundary_u,
			   int* isDiffusiveFluxBoundary_v,
			   int* isDiffusiveFluxBoundary_w,
			   double* ebqe_bc_p_ext,
			   double* ebqe_bc_flux_mass_ext,
			   double* ebqe_bc_flux_mom_u_adv_ext,
			   double* ebqe_bc_flux_mom_v_adv_ext,
			   double* ebqe_bc_flux_mom_w_adv_ext,
			   double* ebqe_bc_u_ext,
			   double* ebqe_bc_flux_u_diff_ext,
			   double* ebqe_penalty_ext,
			   double* ebqe_bc_v_ext,
			   double* ebqe_bc_flux_v_diff_ext,
			   double* ebqe_bc_w_ext,
			   double* ebqe_bc_flux_w_diff_ext,
			   int* csrColumnOffsets_eb_p_p,
			   int* csrColumnOffsets_eb_p_u,
			   int* csrColumnOffsets_eb_p_v,
			   int* csrColumnOffsets_eb_p_w,
			   int* csrColumnOffsets_eb_u_p,
			   int* csrColumnOffsets_eb_u_u,
			   int* csrColumnOffsets_eb_u_v,
			   int* csrColumnOffsets_eb_u_w,
			   int* csrColumnOffsets_eb_v_p,
			   int* csrColumnOffsets_eb_v_u,
			   int* csrColumnOffsets_eb_v_v,
			   int* csrColumnOffsets_eb_v_w,
			   int* csrColumnOffsets_eb_w_p,
			   int* csrColumnOffsets_eb_w_u,
			   int* csrColumnOffsets_eb_w_v,
			   int* csrColumnOffsets_eb_w_w)
    {
      //
      //loop over elements to compute volume integrals and load them into the element Jacobians and global Jacobian
      //
      for(int eN=0;eN<nElements_global;eN++)
	{
	  register double eps_rho,eps_mu;

	  register double  elementJacobian_p_p[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_p_u[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_p_v[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_p_w[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_u_p[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_u_u[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_u_v[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_u_w[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_v_p[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_v_u[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_v_v[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_v_w[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_w_p[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_w_u[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_w_v[nDOF_test_element][nDOF_trial_element],
	    elementJacobian_w_w[nDOF_test_element][nDOF_trial_element];
	  for (int i=0;i<nDOF_test_element;i++)
	    for (int j=0;j<nDOF_trial_element;j++)
	      {
		elementJacobian_p_p[i][j]=0.0;
		elementJacobian_p_u[i][j]=0.0;
		elementJacobian_p_v[i][j]=0.0;
		elementJacobian_p_w[i][j]=0.0;
		elementJacobian_u_p[i][j]=0.0;
		elementJacobian_u_u[i][j]=0.0;
		elementJacobian_u_v[i][j]=0.0;
		elementJacobian_u_w[i][j]=0.0;
		elementJacobian_v_p[i][j]=0.0;
		elementJacobian_v_u[i][j]=0.0;
		elementJacobian_v_v[i][j]=0.0;
		elementJacobian_v_w[i][j]=0.0;
		elementJacobian_w_p[i][j]=0.0;
		elementJacobian_w_u[i][j]=0.0;
		elementJacobian_w_v[i][j]=0.0;
		elementJacobian_w_w[i][j]=0.0;
	      }
	  for  (int k=0;k<nQuadraturePoints_element;k++)
	    {
	      int eN_k = eN*nQuadraturePoints_element+k, //index to a scalar at a quadrature point
		eN_k_nSpace = eN_k*nSpace,
		eN_nDOF_trial_element = eN*nDOF_trial_element; //index to a vector at a quadrature point

	      //declare local storage
	      register double p=0.0,C=0.0,u=0.0,v=0.0,w=0.0,
		grad_p[nSpace],grad_C[nSpace],grad_u[nSpace],grad_v[nSpace],grad_w[nSpace],
		mass_acc=0.0,
		dmass_acc=0.0,
		mom_u_acc=0.0,
		dmom_u_acc_u=0.0,
		mom_v_acc=0.0,
		dmom_v_acc_v=0.0,
		mom_w_acc=0.0,
		dmom_w_acc_w=0.0,
		mass_adv[nSpace],
		dmass_adv_u[nSpace],
		dmass_adv_v[nSpace],
		dmass_adv_w[nSpace],
		mom_u_adv[nSpace],
		dmom_u_adv_u[nSpace],
		dmom_u_adv_v[nSpace],
		dmom_u_adv_w[nSpace],
		mom_v_adv[nSpace],
		dmom_v_adv_u[nSpace],
		dmom_v_adv_v[nSpace],
		dmom_v_adv_w[nSpace],
		mom_w_adv[nSpace],
		dmom_w_adv_u[nSpace],
		dmom_w_adv_v[nSpace],
		dmom_w_adv_w[nSpace],
		mom_u_diff_ten[nSpace],
		mom_v_diff_ten[nSpace],
		mom_w_diff_ten[nSpace],
		mom_uv_diff_ten[1],
		mom_uw_diff_ten[1],
		mom_vu_diff_ten[1],
		mom_vw_diff_ten[1],
		mom_wu_diff_ten[1],
		mom_wv_diff_ten[1],
		mom_u_source=0.0,
		dmom_u_source_p=0.0,
		mom_v_source=0.0,
		dmom_v_source_p=0.0,
		mom_w_source=0.0,
		dmom_w_source_p=0.0,
		mom_u_ham=0.0,
		dmom_u_ham_grad_p[nSpace],
		mom_v_ham=0.0,
		dmom_v_ham_grad_p[nSpace],
		mom_w_ham=0.0,
		dmom_w_ham_grad_p[nSpace],
		mass_acc_t=0.0,
		dmass_acc_t=0.0,
		mom_u_acc_t=0.0,
		dmom_u_acc_u_t=0.0,
		mom_v_acc_t=0.0,
		dmom_v_acc_v_t=0.0,
		mom_w_acc_t=0.0,
		dmom_w_acc_w_t=0.0,
		pdeResidual_p=0.0,
		pdeResidual_u=0.0,
		pdeResidual_v=0.0,
		pdeResidual_w=0.0,	    
		dpdeResidual_p_u[nDOF_trial_element],dpdeResidual_p_v[nDOF_trial_element],dpdeResidual_p_w[nDOF_trial_element],
		dpdeResidual_u_p[nDOF_trial_element],dpdeResidual_u_u[nDOF_trial_element],
		dpdeResidual_v_p[nDOF_trial_element],dpdeResidual_v_v[nDOF_trial_element],
		dpdeResidual_w_p[nDOF_trial_element],dpdeResidual_w_w[nDOF_trial_element],
		Lstar_u_p[nDOF_test_element],
		Lstar_v_p[nDOF_test_element],
		Lstar_w_p[nDOF_test_element],
		Lstar_u_u[nDOF_test_element],
		Lstar_v_v[nDOF_test_element],
		Lstar_w_w[nDOF_test_element],
		Lstar_p_u[nDOF_test_element],
		Lstar_p_v[nDOF_test_element],
		Lstar_p_w[nDOF_test_element],
		subgridError_p=0.0,
		subgridError_u=0.0,
		subgridError_v=0.0,
		subgridError_w=0.0,	    
		dsubgridError_p_u[nDOF_trial_element],
		dsubgridError_p_v[nDOF_trial_element],
		dsubgridError_p_w[nDOF_trial_element],
		dsubgridError_u_p[nDOF_trial_element],
		dsubgridError_u_u[nDOF_trial_element],
		dsubgridError_v_p[nDOF_trial_element],
		dsubgridError_v_v[nDOF_trial_element],
		dsubgridError_w_p[nDOF_trial_element],
		dsubgridError_w_w[nDOF_trial_element],
		tau_p=0.0,tau_p0=0.0,tau_p1=0.0,
		tau_v=0.0,tau_v0=0.0,tau_v1=0.0,
		jac[nSpace*nSpace],
		jacDet,
		jacInv[nSpace*nSpace],
		p_grad_trial[nDOF_trial_element*nSpace],vel_grad_trial[nDOF_trial_element*nSpace],
		dV,
		p_test_dV[nDOF_test_element],vel_test_dV[nDOF_test_element],
		p_grad_test_dV[nDOF_test_element*nSpace],vel_grad_test_dV[nDOF_test_element*nSpace],
		x,y,z,xt,yt,zt,
		G[nSpace*nSpace],G_dd_G,tr_G,h_phi, dmom_adv_star[nSpace], dmom_adv_sge[nSpace];
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
	      //xt=0.0;yt=0.0;zt=0.0;
	      //std::cout<<"xt "<<xt<<'\t'<<yt<<'\t'<<zt<<std::endl;
	      //get the physical integration weight
	      dV = fabs(jacDet)*dV_ref[k];
	      ck.calculateG(jacInv,G,G_dd_G,tr_G);
	      ck.calculateGScale(G,&normal_phi[eN_k_nSpace],h_phi);
	
	      eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      
	      //get the trial function gradients
	      ck.gradTrialFromRef(&p_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,p_grad_trial);
	      ck.gradTrialFromRef(&vel_grad_trial_ref[k*nDOF_trial_element*nSpace],jacInv,vel_grad_trial);
	      //get the solution 	
	      ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],p);
	      ck.valFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_ref[k*nDOF_trial_element],C);
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],u);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],v);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_ref[k*nDOF_trial_element],w);
	      //get the solution gradients
	      ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_p);
	      ck.gradFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial,grad_C);
	      ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_u);
	      ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_v);
	      ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial,grad_w);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  p_test_dV[j] = p_test_ref[k*nDOF_trial_element+j]*dV;
		  vel_test_dV[j] = vel_test_ref[k*nDOF_trial_element+j]*dV;
		  for (int I=0;I<nSpace;I++)
		    {
		      p_grad_test_dV[j*nSpace+I]   = p_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin
		      vel_grad_test_dV[j*nSpace+I] = vel_grad_trial[j*nSpace+I]*dV;//cek warning won't work for Petrov-Galerkin}
		    }
		}
	      //cek debug
	      // q_velocity_sge[eN_k_nSpace+0]=1.0;
	      // q_velocity_sge[eN_k_nSpace+1]=1.0;
	      // q_velocity_sge[eN_k_nSpace+2]=1.0;
	      // //
	      // //debugging section for finite element calculations on interior
	      // //
	      // std::cout<<"eN = "<<eN<<" k = "<<k<<std::endl;
	      // for (int j=0;j<nDOF_trial_element;j++)
	      //   {
	      //     std::cout<<"p_trial["<<j<<"] "<<p_trial_ref[k*nDOF_test_element+j]<<std::endl;
	      //     std::cout<<"vel_trial["<<j<<"] "<<vel_trial_ref[k*nDOF_test_element+j]<<std::endl;
	      //     for (int I=0;I<nSpace;I++)
	      // 	{
	      // 	  std::cout<<"p_grad_trial["<<j<<","<<I<<"] "<<p_grad_trial[j*nSpace+I]<<std::endl;
	      // 	  std::cout<<"vel_grad_trial["<<j<<","<<I<<"] "<<vel_grad_trial[j*nSpace+I]<<std::endl;
	      // 	  std::cout<<"p_grad_test_dV["<<j<<","<<I<<"] "<<p_grad_test_dV[j*nSpace+I]<<std::endl;
	      // 	  std::cout<<"vel_grad_test_dV["<<j<<","<<I<<"] "<<vel_grad_test_dV[j*nSpace+I]<<std::endl;
	      // 	}
	      //   }
	      // std::cout<<"p "<<p<<std::endl;
	      // std::cout<<"u "<<u<<std::endl;
	      // std::cout<<"v "<<v<<std::endl;
	      // std::cout<<"w "<<w<<std::endl;
	      // for(int I=0;I<nSpace;I++)
	      //   {
	      //     std::cout<<"grad_p["<<I<<"] "<<grad_p[I]<<std::endl;
	      //     std::cout<<"grad_u["<<I<<"] "<<grad_u[I]<<std::endl;
	      //     std::cout<<"grad_v["<<I<<"] "<<grad_v[I]<<std::endl;
	      //     std::cout<<"grad_w["<<I<<"] "<<grad_w[I]<<std::endl;
	      //   }
	      //
	      //calculate pde coefficients and derivatives at quadrature points
	      //
	      //cek debug
	      //eps_rho=0.1;
	      //eps_mu=0.1;
	      evaluateCoefficients(eps_rho,
				   eps_mu,
				   sigma,
				   rho_0,
				   nu_0,
				   rho_1,
				   nu_1,
				   g,
				   &q_sedimentForce[eN_k_nSpace],
				   phi[eN_k],
				   &normal_phi[eN_k_nSpace],
				   kappa_phi[eN_k],
				   C,
				   grad_C,
				   p,
				   grad_p,
				   u,
				   v,
				   w,
				   mass_acc,
				   mom_u_acc,
				   dmom_u_acc_u,
				   mom_v_acc,
				   dmom_v_acc_v,
				   mom_w_acc,
				   dmom_w_acc_w,
				   mass_adv,
				   dmass_adv_u,
				   dmass_adv_v,
				   dmass_adv_w,
				   mom_u_adv,
				   dmom_u_adv_u,
				   dmom_u_adv_v,
				   dmom_u_adv_w,
				   mom_v_adv,
				   dmom_v_adv_u,
				   dmom_v_adv_v,
				   dmom_v_adv_w,
				   mom_w_adv,
				   dmom_w_adv_u,
				   dmom_w_adv_v,
				   dmom_w_adv_w,
				   mom_u_diff_ten,
				   mom_v_diff_ten,
				   mom_w_diff_ten,
				   mom_uv_diff_ten,
				   mom_uw_diff_ten,
				   mom_vu_diff_ten,
				   mom_vw_diff_ten,
				   mom_wu_diff_ten,
				   mom_wv_diff_ten,
				   mom_u_source,
				   dmom_u_source_p,
				   mom_v_source,
				   dmom_v_source_p,
				   mom_w_source,
				   dmom_w_source_p,
				   mom_u_ham,
				   dmom_u_ham_grad_p,
				   mom_v_ham,
				   dmom_v_ham_grad_p,
				   mom_w_ham,
				   dmom_w_ham_grad_p);          
	      //
	      //moving mesh
	      //
	      mass_adv[0] -= MOVING_DOMAIN*mass_acc*xt;
	      mass_adv[1] -= MOVING_DOMAIN*mass_acc*yt;
	      mass_adv[2] -= MOVING_DOMAIN*mass_acc*zt;

	      mom_u_adv[0] -= MOVING_DOMAIN*mom_u_acc*xt;
	      mom_u_adv[1] -= MOVING_DOMAIN*mom_u_acc*yt;
	      mom_u_adv[2] -= MOVING_DOMAIN*mom_u_acc*zt;
	      dmom_u_adv_u[0] -= MOVING_DOMAIN*dmom_u_acc_u*xt;
	      dmom_u_adv_u[1] -= MOVING_DOMAIN*dmom_u_acc_u*yt;
	      dmom_u_adv_u[2] -= MOVING_DOMAIN*dmom_u_acc_u*zt;

	      mom_v_adv[0] -= MOVING_DOMAIN*mom_v_acc*xt;
	      mom_v_adv[1] -= MOVING_DOMAIN*mom_v_acc*yt;
	      mom_v_adv[2] -= MOVING_DOMAIN*mom_v_acc*zt;
	      dmom_v_adv_v[0] -= MOVING_DOMAIN*dmom_v_acc_v*xt;
	      dmom_v_adv_v[1] -= MOVING_DOMAIN*dmom_v_acc_v*yt;
	      dmom_v_adv_v[2] -= MOVING_DOMAIN*dmom_v_acc_v*zt;

	      mom_w_adv[0] -= MOVING_DOMAIN*mom_w_acc*xt;
	      mom_w_adv[1] -= MOVING_DOMAIN*mom_w_acc*yt;
	      mom_w_adv[2] -= MOVING_DOMAIN*mom_w_acc*zt;
	      dmom_w_adv_w[0] -= MOVING_DOMAIN*dmom_w_acc_w*xt;
	      dmom_w_adv_w[1] -= MOVING_DOMAIN*dmom_w_acc_w*yt;
	      dmom_w_adv_w[2] -= MOVING_DOMAIN*dmom_w_acc_w*zt;
	      //
	      //calculate time derivatives
	      //
	      ck.bdf(alphaBDF,
		     q_mass_acc_beta_bdf[eN_k],
		     mass_acc,
		     dmass_acc,
		     mass_acc_t,
		     dmass_acc_t);
	      ck.bdf(alphaBDF,
		     q_mom_u_acc_beta_bdf[eN_k],
		     mom_u_acc,
		     dmom_u_acc_u,
		     mom_u_acc_t,
		     dmom_u_acc_u_t);
	      ck.bdf(alphaBDF,
		     q_mom_v_acc_beta_bdf[eN_k],
		     mom_v_acc,
		     dmom_v_acc_v,
		     mom_v_acc_t,
		     dmom_v_acc_v_t);
	      ck.bdf(alphaBDF,
		     q_mom_w_acc_beta_bdf[eN_k],
		     mom_w_acc,
		     dmom_w_acc_w,
		     mom_w_acc_t,
		     dmom_w_acc_w_t);
	      //
	      //calculate subgrid error contribution to the Jacobian (strong residual, adjoint, jacobian of strong residual)
	      //
              dmom_adv_sge[0] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+0] - MOVING_DOMAIN*xt);
              dmom_adv_sge[1] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+1] - MOVING_DOMAIN*yt);
              dmom_adv_sge[2] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+2] - MOVING_DOMAIN*zt);


	      //
	      //calculate strong residual
	      //
	      pdeResidual_p = ck.Mass_strong(mass_acc_t) +
		ck.Advection_strong(dmass_adv_u,grad_u) +
		ck.Advection_strong(dmass_adv_v,grad_v) +
		ck.Advection_strong(dmass_adv_w,grad_w);
	      
	      pdeResidual_u = ck.Mass_strong(mom_u_acc_t) +
		ck.Advection_strong(dmom_adv_sge,grad_u) +
		ck.Hamiltonian_strong(dmom_u_ham_grad_p,grad_p) +
		ck.Reaction_strong(mom_u_source);
	  
	      pdeResidual_v = ck.Mass_strong(mom_v_acc_t) +
		ck.Advection_strong(dmom_adv_sge,grad_v) +
		ck.Hamiltonian_strong(dmom_v_ham_grad_p,grad_p) + 
		ck.Reaction_strong(mom_v_source);
	  
	      pdeResidual_w = ck.Mass_strong(mom_w_acc_t) + 
		ck.Advection_strong(dmom_adv_sge,grad_w) +
		ck.Hamiltonian_strong(dmom_w_ham_grad_p,grad_p) +
		ck.Reaction_strong(mom_w_source);

	      //calculate the Jacobian of strong residual
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  register int j_nSpace = j*nSpace;
		  dpdeResidual_p_u[j]=ck.AdvectionJacobian_strong(dmass_adv_u,&vel_grad_trial[j_nSpace]);
		  dpdeResidual_p_v[j]=ck.AdvectionJacobian_strong(dmass_adv_v,&vel_grad_trial[j_nSpace]);
		  dpdeResidual_p_w[j]=ck.AdvectionJacobian_strong(dmass_adv_w,&vel_grad_trial[j_nSpace]);

		  dpdeResidual_u_p[j]=ck.HamiltonianJacobian_strong(dmom_u_ham_grad_p,&p_grad_trial[j_nSpace]) +
		    ck.ReactionJacobian_strong(dmom_u_source_p,p_trial_ref[k*nDOF_trial_element+j]);
		  dpdeResidual_u_u[j]=ck.MassJacobian_strong(dmom_u_acc_u_t,vel_trial_ref[k*nDOF_trial_element+j]) +
		    ck.AdvectionJacobian_strong(dmom_adv_sge,&vel_grad_trial[j_nSpace]);
	      
		  dpdeResidual_v_p[j]=ck.HamiltonianJacobian_strong(dmom_v_ham_grad_p,&p_grad_trial[j_nSpace]) +
		    ck.ReactionJacobian_strong(dmom_v_source_p,p_trial_ref[k*nDOF_trial_element+j]);
		  dpdeResidual_v_v[j]=ck.MassJacobian_strong(dmom_v_acc_v_t,vel_trial_ref[k*nDOF_trial_element+j]) +
		    ck.AdvectionJacobian_strong(dmom_adv_sge,&vel_grad_trial[j_nSpace]);
	      
		  dpdeResidual_w_p[j]=ck.HamiltonianJacobian_strong(dmom_w_ham_grad_p,&p_grad_trial[j_nSpace]) +
		    ck.ReactionJacobian_strong(dmom_w_source_p,p_trial_ref[k*nDOF_trial_element+j]);
		  dpdeResidual_w_w[j]=ck.MassJacobian_strong(dmom_w_acc_w_t,vel_trial_ref[k*nDOF_trial_element+j]) + 
		    ck.AdvectionJacobian_strong(dmom_adv_sge,&vel_grad_trial[j_nSpace]);
		}
	      //calculate tau and tau*Res
	      calculateSubgridError_tau(hFactor,
					elementDiameter[eN],
					dmom_u_acc_u_t,
					dmom_u_acc_u,
					dmom_adv_sge,
					mom_u_diff_ten[1],
					dmom_u_ham_grad_p[0],
					tau_v0,
					tau_p0,
					q_cfl[eN_k]);
					
	      calculateSubgridError_tau(Ct_sge,Cd_sge,
			                G,G_dd_G,tr_G,
					dmom_u_acc_u_t,
					dmom_adv_sge,
					mom_u_diff_ten[1],
                                        dmom_u_ham_grad_p[0],					
					tau_v1,
					tau_p1,
					q_cfl[eN_k]);					
					
					
	      tau_v = useMetrics*tau_v1+(1.0-useMetrics)*tau_v0;
	      tau_p = useMetrics*tau_p1+(1.0-useMetrics)*tau_p0;					
					
					
	      //calculateSubgridError_tau(Ct_sge,
	      // 				      Cd_sge,
	      // 				      G,
	      // 				      G_dd_G,
	      // 				      tr_G,
	      // 				      dmom_u_acc_u,//rho
	      // 				      dmom_u_acc_u_t/dmom_u_acc_u,//Dt
	      // 				      &q_velocity_sge[eN_k_nSpace],//v
	      // 				      mom_u_diff_ten[1],//mu
	      // 				      tau_p,
	      // 				      tau_v,
	      // 				      q_cfl[eN_k]);
	      //cek debug
	      //tau_p = 0.0;
	      //tau_v = 0.0;
	      
	      
	      calculateSubgridError_tauRes(tau_p,
					   tau_v,
					   pdeResidual_p,
					   pdeResidual_u,
					   pdeResidual_v,
					   pdeResidual_w,
					   subgridError_p,
					   subgridError_u,
					   subgridError_v,
					   subgridError_w);	      
	      
	      calculateSubgridErrorDerivatives_tauRes(tau_p,
						      tau_v,
						      dpdeResidual_p_u,
						      dpdeResidual_p_v,
						      dpdeResidual_p_w,
						      dpdeResidual_u_p,
						      dpdeResidual_u_u,
						      dpdeResidual_v_p,
						      dpdeResidual_v_v,
						      dpdeResidual_w_p,
						      dpdeResidual_w_w,
						      dsubgridError_p_u,
						      dsubgridError_p_v,
						      dsubgridError_p_w,
						      dsubgridError_u_p,
						      dsubgridError_u_u,
						      dsubgridError_v_p,
						      dsubgridError_v_v,
						      dsubgridError_w_p,
						      dsubgridError_w_w);
	      // velocity used in adjoint (VMS or RBLES, with or without lagging the grid scale velocity)
	      dmom_adv_star[0] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+0] - MOVING_DOMAIN*xt + useRBLES*subgridError_u);
	      dmom_adv_star[1] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+1] - MOVING_DOMAIN*yt + useRBLES*subgridError_v);
	      dmom_adv_star[2] = dmom_u_acc_u*(q_velocity_sge[eN_k_nSpace+2] - MOVING_DOMAIN*zt + useRBLES*subgridError_w);
          
	      //calculate the adjoint times the test functions
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int i_nSpace = i*nSpace;
		  Lstar_u_p[i]=ck.Advection_adjoint(dmass_adv_u,&p_grad_test_dV[i_nSpace]);
		  Lstar_v_p[i]=ck.Advection_adjoint(dmass_adv_v,&p_grad_test_dV[i_nSpace]);
		  Lstar_w_p[i]=ck.Advection_adjoint(dmass_adv_w,&p_grad_test_dV[i_nSpace]);
		  Lstar_u_u[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
		  Lstar_v_v[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
		  Lstar_w_w[i]=ck.Advection_adjoint(dmom_adv_star,&vel_grad_test_dV[i_nSpace]);
		  Lstar_p_u[i]=ck.Hamiltonian_adjoint(dmom_u_ham_grad_p,&vel_grad_test_dV[i_nSpace]) + ck.Reaction_adjoint(dmom_u_source_p,vel_test_dV[i]);
		  Lstar_p_v[i]=ck.Hamiltonian_adjoint(dmom_v_ham_grad_p,&vel_grad_test_dV[i_nSpace]) + ck.Reaction_adjoint(dmom_v_source_p,vel_test_dV[i]);
		  Lstar_p_w[i]=ck.Hamiltonian_adjoint(dmom_w_ham_grad_p,&vel_grad_test_dV[i_nSpace]) + ck.Reaction_adjoint(dmom_w_source_p,vel_test_dV[i]);
		}


	      //cek todo add RBLES terms consistent to residual modifications or ignore them partials w.r.t the additional RBLES terms
	      for(int i=0;i<nDOF_test_element;i++)
		{
		  register int i_nSpace = i*nSpace;
		  for(int j=0;j<nDOF_trial_element;j++) 
		    { 
		      register int j_nSpace = j*nSpace;
		      elementJacobian_p_p[i][j] += ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_p[i]) + 
			ck.SubgridErrorJacobian(dsubgridError_v_p[j],Lstar_v_p[i]) + 
			ck.SubgridErrorJacobian(dsubgridError_w_p[j],Lstar_w_p[i]); 

		      elementJacobian_p_u[i][j] += ck.AdvectionJacobian_weak(dmass_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u_p[i]); 
		      elementJacobian_p_v[i][j] += ck.AdvectionJacobian_weak(dmass_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_v_v[j],Lstar_v_p[i]); 
		      elementJacobian_p_w[i][j] += ck.AdvectionJacobian_weak(dmass_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&p_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_w_w[j],Lstar_w_p[i]); 

		      elementJacobian_u_p[i][j] += ck.HamiltonianJacobian_weak(dmom_u_ham_grad_p,&p_grad_trial[j_nSpace],vel_test_dV[i]) + 
			ck.ReactionJacobian_weak(dmom_u_source_p,p_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.SubgridErrorJacobian(dsubgridError_u_p[j],Lstar_u_u[i]); 
		      elementJacobian_u_u[i][j] += ck.MassJacobian_weak(dmom_u_acc_u_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.AdvectionJacobian_weak(dmom_u_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			ck.SimpleDiffusionJacobian_weak(sdInfo_u_u_rowptr,sdInfo_u_u_colind,mom_u_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_u[i]) + 
			ck.SubgridErrorJacobian(dsubgridError_u_u[j],Lstar_u_u[i]) + 
			ck.NumericalDiffusionJacobian(q_numDiff_u_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]); 
		      elementJacobian_u_v[i][j] += ck.AdvectionJacobian_weak(dmom_u_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
			ck.SimpleDiffusionJacobian_weak(sdInfo_u_v_rowptr,sdInfo_u_v_colind,mom_uv_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_u[i]); 
		      elementJacobian_u_w[i][j] += ck.AdvectionJacobian_weak(dmom_u_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
			ck.SimpleDiffusionJacobian_weak(sdInfo_u_w_rowptr,sdInfo_u_w_colind,mom_uw_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_w[j],Lstar_p_u[i]); 

		      elementJacobian_v_p[i][j] += ck.HamiltonianJacobian_weak(dmom_v_ham_grad_p,&p_grad_trial[j_nSpace],vel_test_dV[i]) + 
			ck.ReactionJacobian_weak(dmom_v_source_p,p_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.SubgridErrorJacobian(dsubgridError_v_p[j],Lstar_v_v[i]); 
		      elementJacobian_v_u[i][j] += ck.AdvectionJacobian_weak(dmom_v_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
			ck.SimpleDiffusionJacobian_weak(sdInfo_v_u_rowptr,sdInfo_v_u_colind,mom_vu_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_v[i]);
		      elementJacobian_v_v[i][j] += ck.MassJacobian_weak(dmom_v_acc_v_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.AdvectionJacobian_weak(dmom_v_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +
			ck.SimpleDiffusionJacobian_weak(sdInfo_v_v_rowptr,sdInfo_v_v_colind,mom_v_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_v[i]) +
			ck.SubgridErrorJacobian(dsubgridError_v_v[j],Lstar_v_v[i]) + 
			ck.NumericalDiffusionJacobian(q_numDiff_v_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]); 
		      elementJacobian_v_w[i][j] += ck.AdvectionJacobian_weak(dmom_v_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +  
			ck.SimpleDiffusionJacobian_weak(sdInfo_v_w_rowptr,sdInfo_v_w_colind,mom_vw_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_w[j],Lstar_p_v[i]);

		      elementJacobian_w_p[i][j] += ck.HamiltonianJacobian_weak(dmom_w_ham_grad_p,&p_grad_trial[j_nSpace],vel_test_dV[i]) + 
			ck.ReactionJacobian_weak(dmom_w_source_p,p_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.SubgridErrorJacobian(dsubgridError_w_p[j],Lstar_w_w[i]); 
		      elementJacobian_w_u[i][j] += ck.AdvectionJacobian_weak(dmom_w_adv_u,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +  
			ck.SimpleDiffusionJacobian_weak(sdInfo_w_u_rowptr,sdInfo_w_u_colind,mom_wu_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_u[j],Lstar_p_w[i]); 
		      elementJacobian_w_v[i][j] += ck.AdvectionJacobian_weak(dmom_w_adv_v,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) + 
			ck.SimpleDiffusionJacobian_weak(sdInfo_w_v_rowptr,sdInfo_w_v_colind,mom_wv_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_v[j],Lstar_p_w[i]); 
		      elementJacobian_w_w[i][j] += ck.MassJacobian_weak(dmom_w_acc_w_t,vel_trial_ref[k*nDOF_trial_element+j],vel_test_dV[i]) + 
			ck.AdvectionJacobian_weak(dmom_w_adv_w,vel_trial_ref[k*nDOF_trial_element+j],&vel_grad_test_dV[i_nSpace]) +  
			ck.SimpleDiffusionJacobian_weak(sdInfo_w_w_rowptr,sdInfo_w_w_colind,mom_w_diff_ten,&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]) + 
			ck.SubgridErrorJacobian(dsubgridError_p_w[j],Lstar_p_w[i]) + 
			ck.SubgridErrorJacobian(dsubgridError_w_w[j],Lstar_w_w[i]) + 
			ck.NumericalDiffusionJacobian(q_numDiff_w_last[eN_k],&vel_grad_trial[j_nSpace],&vel_grad_test_dV[i_nSpace]); 
		    }//j
		}//i
	    }//k
	  //
	  //load into element Jacobian into global Jacobian
	  //
	  for (int i=0;i<nDOF_test_element;i++)
	    {
	      register int eN_i = eN*nDOF_test_element+i;
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  register int eN_i_j = eN_i*nDOF_trial_element+j;
		  globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_p_p[eN_i_j]] += elementJacobian_p_p[i][j];
		  globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_p_u[eN_i_j]] += elementJacobian_p_u[i][j];
		  globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_p_v[eN_i_j]] += elementJacobian_p_v[i][j];
		  globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_p_w[eN_i_j]] += elementJacobian_p_w[i][j];

		  globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_u_p[eN_i_j]] += elementJacobian_u_p[i][j];
		  globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_u_u[eN_i_j]] += elementJacobian_u_u[i][j];
		  globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_u_v[eN_i_j]] += elementJacobian_u_v[i][j];
		  globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_u_w[eN_i_j]] += elementJacobian_u_w[i][j];

		  globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_v_p[eN_i_j]] += elementJacobian_v_p[i][j];
		  globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_v_u[eN_i_j]] += elementJacobian_v_u[i][j];
		  globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_v_v[eN_i_j]] += elementJacobian_v_v[i][j];
		  globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_v_w[eN_i_j]] += elementJacobian_v_w[i][j];

		  globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_w_p[eN_i_j]] += elementJacobian_w_p[i][j];
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
	  register int ebN = exteriorElementBoundariesArray[ebNE],
	    eN  = elementBoundaryElementsArray[ebN*2+0],
	    eN_nDOF_trial_element = eN*nDOF_trial_element,
	    ebN_local = elementBoundaryLocalElementBoundariesArray[ebN*2+0];
	  register double eps_rho,eps_mu;
	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++) 
	    { 
	      register int ebNE_kb = ebNE*nQuadraturePoints_elementBoundary+kb,
		ebNE_kb_nSpace = ebNE_kb*nSpace,
		ebN_local_kb = ebN_local*nQuadraturePoints_elementBoundary+kb,
		ebN_local_kb_nSpace = ebN_local_kb*nSpace;

	      register double p_ext=0.0,C_ext=0.0,
		u_ext=0.0,
		v_ext=0.0,
		w_ext=0.0,
		grad_p_ext[nSpace],
		grad_C_ext[nSpace],
		grad_u_ext[nSpace],
		grad_v_ext[nSpace],
		grad_w_ext[nSpace],
		mass_acc_ext=0.0,
		dmass_acc_ext=0.0,
		mom_u_acc_ext=0.0,
		dmom_u_acc_u_ext=0.0,
		mom_v_acc_ext=0.0,
		dmom_v_acc_v_ext=0.0,
		mom_w_acc_ext=0.0,
		dmom_w_acc_w_ext=0.0,
		mass_adv_ext[nSpace],
		dmass_adv_u_ext[nSpace],
		dmass_adv_v_ext[nSpace],
		dmass_adv_w_ext[nSpace],
		mom_u_adv_ext[nSpace],
		dmom_u_adv_u_ext[nSpace],
		dmom_u_adv_v_ext[nSpace],
		dmom_u_adv_w_ext[nSpace],
		mom_v_adv_ext[nSpace],
		dmom_v_adv_u_ext[nSpace],
		dmom_v_adv_v_ext[nSpace],
		dmom_v_adv_w_ext[nSpace],
		mom_w_adv_ext[nSpace],
		dmom_w_adv_u_ext[nSpace],
		dmom_w_adv_v_ext[nSpace],
		dmom_w_adv_w_ext[nSpace],
		mom_u_diff_ten_ext[nSpace],
		mom_v_diff_ten_ext[nSpace],
		mom_w_diff_ten_ext[nSpace],
		mom_uv_diff_ten_ext[1],
		mom_uw_diff_ten_ext[1],
		mom_vu_diff_ten_ext[1],
		mom_vw_diff_ten_ext[1],
		mom_wu_diff_ten_ext[1],
		mom_wv_diff_ten_ext[1],
		mom_u_source_ext=0.0,
		dmom_u_source_p_ext=0.0,
		mom_v_source_ext=0.0,
		dmom_v_source_p_ext=0.0,
		mom_w_source_ext=0.0,
		dmom_w_source_p_ext=0.0,
		mom_u_ham_ext=0.0,
		dmom_u_ham_grad_p_ext[nSpace],
		mom_v_ham_ext=0.0,
		dmom_v_ham_grad_p_ext[nSpace],
		mom_w_ham_ext=0.0,
		dmom_w_ham_grad_p_ext[nSpace],
		dmom_u_adv_p_ext[nSpace],
		dmom_v_adv_p_ext[nSpace],
		dmom_w_adv_p_ext[nSpace],
		dflux_mass_u_ext=0.0,
		dflux_mass_v_ext=0.0,
		dflux_mass_w_ext=0.0,
		dflux_mom_u_adv_p_ext=0.0,
		dflux_mom_u_adv_u_ext=0.0,
		dflux_mom_u_adv_v_ext=0.0,
		dflux_mom_u_adv_w_ext=0.0,
		dflux_mom_v_adv_p_ext=0.0,
		dflux_mom_v_adv_u_ext=0.0,
		dflux_mom_v_adv_v_ext=0.0,
		dflux_mom_v_adv_w_ext=0.0,
		dflux_mom_w_adv_p_ext=0.0,
		dflux_mom_w_adv_u_ext=0.0,
		dflux_mom_w_adv_v_ext=0.0,
		dflux_mom_w_adv_w_ext=0.0,
		bc_p_ext=0.0,
		bc_C_ext=0.0,
		bc_u_ext=0.0,
		bc_v_ext=0.0,
		bc_w_ext=0.0,
		bc_mass_acc_ext=0.0,
		bc_dmass_acc_ext=0.0,
		bc_mom_u_acc_ext=0.0,
		bc_dmom_u_acc_u_ext=0.0,
		bc_mom_v_acc_ext=0.0,
		bc_dmom_v_acc_v_ext=0.0,
		bc_mom_w_acc_ext=0.0,
		bc_dmom_w_acc_w_ext=0.0,
		bc_mass_adv_ext[nSpace],
		bc_dmass_adv_u_ext[nSpace],
		bc_dmass_adv_v_ext[nSpace],
		bc_dmass_adv_w_ext[nSpace],
		bc_mom_u_adv_ext[nSpace],
		bc_dmom_u_adv_u_ext[nSpace],
		bc_dmom_u_adv_v_ext[nSpace],
		bc_dmom_u_adv_w_ext[nSpace],
		bc_mom_v_adv_ext[nSpace],
		bc_dmom_v_adv_u_ext[nSpace],
		bc_dmom_v_adv_v_ext[nSpace],
		bc_dmom_v_adv_w_ext[nSpace],
		bc_mom_w_adv_ext[nSpace],
		bc_dmom_w_adv_u_ext[nSpace],
		bc_dmom_w_adv_v_ext[nSpace],
		bc_dmom_w_adv_w_ext[nSpace],
		bc_mom_u_diff_ten_ext[nSpace],
		bc_mom_v_diff_ten_ext[nSpace],
		bc_mom_w_diff_ten_ext[nSpace],
		bc_mom_uv_diff_ten_ext[1],
		bc_mom_uw_diff_ten_ext[1],
		bc_mom_vu_diff_ten_ext[1],
		bc_mom_vw_diff_ten_ext[1],
		bc_mom_wu_diff_ten_ext[1],
		bc_mom_wv_diff_ten_ext[1],
		bc_mom_u_source_ext=0.0,
		bc_mom_v_source_ext=0.0,
		bc_mom_w_source_ext=0.0,
		bc_mom_u_ham_ext=0.0,
		bc_dmom_u_ham_grad_p_ext[nSpace],
		bc_mom_v_ham_ext=0.0,
		bc_dmom_v_ham_grad_p_ext[nSpace],
		bc_mom_w_ham_ext=0.0,
		bc_dmom_w_ham_grad_p_ext[nSpace],
		fluxJacobian_p_p[nDOF_trial_element],
		fluxJacobian_p_u[nDOF_trial_element],
		fluxJacobian_p_v[nDOF_trial_element],
		fluxJacobian_p_w[nDOF_trial_element],
		fluxJacobian_u_p[nDOF_trial_element],
		fluxJacobian_u_u[nDOF_trial_element],
		fluxJacobian_u_v[nDOF_trial_element],
		fluxJacobian_u_w[nDOF_trial_element],
		fluxJacobian_v_p[nDOF_trial_element],
		fluxJacobian_v_u[nDOF_trial_element],
		fluxJacobian_v_v[nDOF_trial_element],
		fluxJacobian_v_w[nDOF_trial_element],
		fluxJacobian_w_p[nDOF_trial_element],
		fluxJacobian_w_u[nDOF_trial_element],
		fluxJacobian_w_v[nDOF_trial_element],
		fluxJacobian_w_w[nDOF_trial_element],
		jac_ext[nSpace*nSpace],
		jacDet_ext,
		jacInv_ext[nSpace*nSpace],
		boundaryJac[nSpace*(nSpace-1)],
		metricTensor[(nSpace-1)*(nSpace-1)],
		metricTensorDetSqrt,
		p_grad_trial_trace[nDOF_trial_element*nSpace],
		vel_grad_trial_trace[nDOF_trial_element*nSpace],
		dS,
		p_test_dS[nDOF_test_element],
		vel_test_dS[nDOF_test_element],
		normal[3],
		x_ext,y_ext,z_ext,xt_ext,yt_ext,zt_ext,integralScaling,
		G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty;
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
	      //xt_ext=0.0;yt_ext=0.0;zt_ext=0.0;
	      //std::cout<<"xt_ext "<<xt_ext<<'\t'<<yt_ext<<'\t'<<zt_ext<<std::endl;
	      dS = ((1.0-MOVING_DOMAIN)*metricTensorDetSqrt + MOVING_DOMAIN*integralScaling)*dS_ref[kb];
	      ck.calculateG(jacInv_ext,G,G_dd_G,tr_G);
	      ck.calculateGScale(G,&ebqe_normal_phi_ext[ebNE_kb_nSpace],h_phi);

	      eps_rho = epsFact_rho*(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);
	      eps_mu  = epsFact_mu *(useMetrics*h_phi+(1.0-useMetrics)*elementDiameter[eN]);

	      //compute shape and solution information
	      //shape
	      ck.gradTrialFromRef(&p_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,p_grad_trial_trace);
	      ck.gradTrialFromRef(&vel_grad_trial_trace_ref[ebN_local_kb_nSpace*nDOF_trial_element],jacInv_ext,vel_grad_trial_trace);
	      //solution and gradients	
	      ck.valFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],p_ext);
	      ck.valFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],&p_trial_trace_ref[ebN_local_kb*nDOF_test_element],C_ext);
	      ck.valFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],u_ext);
	      ck.valFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],v_ext);
	      ck.valFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],&vel_trial_trace_ref[ebN_local_kb*nDOF_test_element],w_ext);
	      ck.gradFromDOF(p_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_p_ext);
	      ck.gradFromDOF(C_dof,&p_l2g[eN_nDOF_trial_element],p_grad_trial_trace,grad_C_ext);
	      ck.gradFromDOF(u_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_u_ext);
	      ck.gradFromDOF(v_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_v_ext);
	      ck.gradFromDOF(w_dof,&vel_l2g[eN_nDOF_trial_element],vel_grad_trial_trace,grad_w_ext);
	      //precalculate test function products with integration weights
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  p_test_dS[j] = p_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		  vel_test_dS[j] = vel_test_trace_ref[ebN_local_kb*nDOF_test_element+j]*dS;
		}
	      // //
	      // //debugging section for finite element calculations on exterior
	      // //
	      // std::cout<<"ebNE = "<<ebNE<<" kb = "<<kb<<std::endl;
	      // for (int j=0;j<nDOF_trial_element;j++)
	      //   {
	      //     std::cout<<"p_trial_trace["<<j<<"] "<<p_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	      //     std::cout<<"vel_trial_trace["<<j<<"] "<<vel_trial_trace_ref[ebN_local_kb*nDOF_trial_element+j]<<std::endl;
	      //     std::cout<<"p_test_dS["<<j<<"] "<<p_test_dS[j]<<std::endl;
	      //     std::cout<<"vel_test_dS["<<j<<"] "<<vel_test_dS[j]<<std::endl;
	      //     for (int I=0;I<nSpace;I++)
	      // 	{
	      // 	  std::cout<<"p_grad_trial_trace["<<j<<","<<I<<"] "<<p_grad_trial_trace[j*nSpace+I]<<std::endl;
	      // 	  std::cout<<"vel_grad_trial_trace["<<j<<","<<I<<"] "<<vel_grad_trial_trace[j*nSpace+I]<<std::endl;
	      // 	}
	      //   }
	      // std::cout<<"p_ext "<<p_ext<<std::endl;
	      // std::cout<<"u_ext "<<u_ext<<std::endl;
	      // std::cout<<"v_ext "<<v_ext<<std::endl;
	      // std::cout<<"w_ext "<<w_ext<<std::endl;
	      // for(int I=0;I<nSpace;I++)
	      //   {
	      //     std::cout<<"grad_p_ext["<<I<<"] "<<grad_p_ext[I]<<std::endl;
	      //     std::cout<<"grad_u_ext["<<I<<"] "<<grad_u_ext[I]<<std::endl;
	      //     std::cout<<"grad_v_ext["<<I<<"] "<<grad_v_ext[I]<<std::endl;
	      //     std::cout<<"grad_w_ext["<<I<<"] "<<grad_w_ext[I]<<std::endl;
	      //   }
	      //
	      //load the boundary values
	      //
	      bc_p_ext = isDOFBoundary_p[ebNE_kb]*ebqe_bc_p_ext[ebNE_kb]+(1-isDOFBoundary_p[ebNE_kb])*p_ext;
	      bc_u_ext = isDOFBoundary_u[ebNE_kb]*ebqe_bc_u_ext[ebNE_kb]+(1-isDOFBoundary_u[ebNE_kb])*u_ext;
	      bc_v_ext = isDOFBoundary_v[ebNE_kb]*ebqe_bc_v_ext[ebNE_kb]+(1-isDOFBoundary_v[ebNE_kb])*v_ext;
	      bc_w_ext = isDOFBoundary_w[ebNE_kb]*ebqe_bc_w_ext[ebNE_kb]+(1-isDOFBoundary_w[ebNE_kb])*w_ext;
	      // 
	      //calculate the internal and external trace of the pde coefficients 
	      // 
	      //cek debug
	      //eps_rho=0.1;
	      //eps_mu=0.1;
	      double ebqe_sedimentForce[nSpace];
	      evaluateCoefficients(eps_rho,
				   eps_mu,
				   sigma,
				   rho_0,
				   nu_0,
				   rho_1,
				   nu_1,
				   g,
				   ebqe_sedimentForce,//[ebNE_kb_nSpace],
				   ebqe_phi_ext[ebNE_kb],
				   &ebqe_normal_phi_ext[ebNE_kb_nSpace],
				   ebqe_kappa_phi_ext[ebNE_kb],
				   C_ext,
				   grad_C_ext,
				   p_ext,
				   grad_p_ext,
				   u_ext,
				   v_ext,
				   w_ext,
				   mass_acc_ext,
				   mom_u_acc_ext,
				   dmom_u_acc_u_ext,
				   mom_v_acc_ext,
				   dmom_v_acc_v_ext,
				   mom_w_acc_ext,
				   dmom_w_acc_w_ext,
				   mass_adv_ext,
				   dmass_adv_u_ext,
				   dmass_adv_v_ext,
				   dmass_adv_w_ext,
				   mom_u_adv_ext,
				   dmom_u_adv_u_ext,
				   dmom_u_adv_v_ext,
				   dmom_u_adv_w_ext,
				   mom_v_adv_ext,
				   dmom_v_adv_u_ext,
				   dmom_v_adv_v_ext,
				   dmom_v_adv_w_ext,
				   mom_w_adv_ext,
				   dmom_w_adv_u_ext,
				   dmom_w_adv_v_ext,
				   dmom_w_adv_w_ext,
				   mom_u_diff_ten_ext,
				   mom_v_diff_ten_ext,
				   mom_w_diff_ten_ext,
				   mom_uv_diff_ten_ext,
				   mom_uw_diff_ten_ext,
				   mom_vu_diff_ten_ext,
				   mom_vw_diff_ten_ext,
				   mom_wu_diff_ten_ext,
				   mom_wv_diff_ten_ext,
				   mom_u_source_ext,
				   dmom_u_source_p_ext,
				   mom_v_source_ext,
				   dmom_v_source_p_ext,
				   mom_w_source_ext,
				   dmom_w_source_p_ext,
				   mom_u_ham_ext,
				   dmom_u_ham_grad_p_ext,
				   mom_v_ham_ext,
				   dmom_v_ham_grad_p_ext,
				   mom_w_ham_ext,
				   dmom_w_ham_grad_p_ext);          
	      evaluateCoefficients(eps_rho,
				   eps_mu,
				   sigma,
				   rho_0,
				   nu_0,
				   rho_1,
				   nu_1,
				   g,
				   ebqe_sedimentForce,//[ebNE_kb_nSpace],
				   ebqe_phi_ext[ebNE_kb],
				   &ebqe_normal_phi_ext[ebNE_kb_nSpace],
				   ebqe_kappa_phi_ext[ebNE_kb],
				   bc_C_ext,
				   grad_C_ext, //cek shouldn't be used
				   bc_p_ext,
				   grad_p_ext, //cek shouldn't be used
				   bc_u_ext,
				   bc_v_ext,
				   bc_w_ext,
				   bc_mass_acc_ext,
				   bc_mom_u_acc_ext,
				   bc_dmom_u_acc_u_ext,
				   bc_mom_v_acc_ext,
				   bc_dmom_v_acc_v_ext,
				   bc_mom_w_acc_ext,
				   bc_dmom_w_acc_w_ext,
				   bc_mass_adv_ext,
				   bc_dmass_adv_u_ext,
				   bc_dmass_adv_v_ext,
				   bc_dmass_adv_w_ext,
				   bc_mom_u_adv_ext,
				   bc_dmom_u_adv_u_ext,
				   bc_dmom_u_adv_v_ext,
				   bc_dmom_u_adv_w_ext,
				   bc_mom_v_adv_ext,
				   bc_dmom_v_adv_u_ext,
				   bc_dmom_v_adv_v_ext,
				   bc_dmom_v_adv_w_ext,
				   bc_mom_w_adv_ext,
				   bc_dmom_w_adv_u_ext,
				   bc_dmom_w_adv_v_ext,
				   bc_dmom_w_adv_w_ext,
				   bc_mom_u_diff_ten_ext,
				   bc_mom_v_diff_ten_ext,
				   bc_mom_w_diff_ten_ext,
				   bc_mom_uv_diff_ten_ext,
				   bc_mom_uw_diff_ten_ext,
				   bc_mom_vu_diff_ten_ext,
				   bc_mom_vw_diff_ten_ext,
				   bc_mom_wu_diff_ten_ext,
				   bc_mom_wv_diff_ten_ext,
				   bc_mom_u_source_ext,
				   dmom_u_source_p_ext,//cek not used
				   bc_mom_v_source_ext,
				   dmom_v_source_p_ext,//cek not used
				   bc_mom_w_source_ext,
				   dmom_w_source_p_ext,//cek not used
				   bc_mom_u_ham_ext,
				   bc_dmom_u_ham_grad_p_ext,
				   bc_mom_v_ham_ext,
				   bc_dmom_v_ham_grad_p_ext,
				   bc_mom_w_ham_ext,
				   bc_dmom_w_ham_grad_p_ext);          
	      //
	      //moving domain
	      //
	      mass_adv_ext[0] -= MOVING_DOMAIN*mass_acc_ext*xt_ext;
	      mass_adv_ext[1] -= MOVING_DOMAIN*mass_acc_ext*yt_ext;
	      mass_adv_ext[2] -= MOVING_DOMAIN*mass_acc_ext*zt_ext;

	      mom_u_adv_ext[0] -= MOVING_DOMAIN*mom_u_acc_ext*xt_ext;
	      mom_u_adv_ext[1] -= MOVING_DOMAIN*mom_u_acc_ext*yt_ext;
	      mom_u_adv_ext[2] -= MOVING_DOMAIN*mom_u_acc_ext*zt_ext;
	      dmom_u_adv_u_ext[0] -= MOVING_DOMAIN*dmom_u_acc_u_ext*xt_ext;
	      dmom_u_adv_u_ext[1] -= MOVING_DOMAIN*dmom_u_acc_u_ext*yt_ext;
	      dmom_u_adv_u_ext[2] -= MOVING_DOMAIN*dmom_u_acc_u_ext*zt_ext;
	      
	      mom_v_adv_ext[0] -= MOVING_DOMAIN*mom_v_acc_ext*xt_ext;
	      mom_v_adv_ext[1] -= MOVING_DOMAIN*mom_v_acc_ext*yt_ext;
	      mom_v_adv_ext[2] -= MOVING_DOMAIN*mom_v_acc_ext*zt_ext;
	      dmom_v_adv_v_ext[0] -= MOVING_DOMAIN*dmom_v_acc_v_ext*xt_ext;
	      dmom_v_adv_v_ext[1] -= MOVING_DOMAIN*dmom_v_acc_v_ext*yt_ext;
	      dmom_v_adv_v_ext[2] -= MOVING_DOMAIN*dmom_v_acc_v_ext*zt_ext;
	      
	      mom_w_adv_ext[0] -= MOVING_DOMAIN*mom_w_acc_ext*xt_ext;
	      mom_w_adv_ext[1] -= MOVING_DOMAIN*mom_w_acc_ext*yt_ext;
	      mom_w_adv_ext[2] -= MOVING_DOMAIN*mom_w_acc_ext*zt_ext;
	      dmom_w_adv_w_ext[0] -= MOVING_DOMAIN*dmom_w_acc_w_ext*xt_ext;
	      dmom_w_adv_w_ext[1] -= MOVING_DOMAIN*dmom_w_acc_w_ext*yt_ext;
	      dmom_w_adv_w_ext[2] -= MOVING_DOMAIN*dmom_w_acc_w_ext*zt_ext;
	      
	      //moving domain bc's
	      bc_mom_u_adv_ext[0] -= MOVING_DOMAIN*bc_mom_u_acc_ext*xt_ext;
	      bc_mom_u_adv_ext[1] -= MOVING_DOMAIN*bc_mom_u_acc_ext*yt_ext;
	      bc_mom_u_adv_ext[2] -= MOVING_DOMAIN*bc_mom_u_acc_ext*zt_ext;
	      
	      bc_mom_v_adv_ext[0] -= MOVING_DOMAIN*bc_mom_v_acc_ext*xt_ext;
	      bc_mom_v_adv_ext[1] -= MOVING_DOMAIN*bc_mom_v_acc_ext*yt_ext;
	      bc_mom_v_adv_ext[2] -= MOVING_DOMAIN*bc_mom_v_acc_ext*zt_ext;

	      bc_mom_w_adv_ext[0] -= MOVING_DOMAIN*bc_mom_w_acc_ext*xt_ext;
	      bc_mom_w_adv_ext[1] -= MOVING_DOMAIN*bc_mom_w_acc_ext*yt_ext;
	      bc_mom_w_adv_ext[2] -= MOVING_DOMAIN*bc_mom_w_acc_ext*zt_ext;
	      // 
	      //calculate the numerical fluxes 
	      // 
	      exteriorNumericalAdvectiveFluxDerivatives(isDOFBoundary_p[ebNE_kb],
							isDOFBoundary_u[ebNE_kb],
							isDOFBoundary_v[ebNE_kb],
							isDOFBoundary_w[ebNE_kb],
							isAdvectiveFluxBoundary_p[ebNE_kb],
							isAdvectiveFluxBoundary_u[ebNE_kb],
							isAdvectiveFluxBoundary_v[ebNE_kb],
							isAdvectiveFluxBoundary_w[ebNE_kb],
							dmom_u_ham_grad_p_ext[0],//=1/rho
							normal,
							bc_p_ext,
							bc_mass_adv_ext,
							bc_mom_u_adv_ext,
							bc_mom_v_adv_ext,
							bc_mom_w_adv_ext,
							ebqe_bc_flux_mass_ext[ebNE_kb],
							ebqe_bc_flux_mom_u_adv_ext[ebNE_kb],
							ebqe_bc_flux_mom_v_adv_ext[ebNE_kb],
							ebqe_bc_flux_mom_w_adv_ext[ebNE_kb],
							p_ext,
							mass_adv_ext,
							mom_u_adv_ext,
							mom_v_adv_ext,
							mom_w_adv_ext,
							dmass_adv_u_ext,
							dmass_adv_v_ext,
							dmass_adv_w_ext,
							dmom_u_adv_p_ext,
							dmom_u_adv_u_ext,
							dmom_u_adv_v_ext,
							dmom_u_adv_w_ext,
							dmom_v_adv_p_ext,
							dmom_v_adv_u_ext,
							dmom_v_adv_v_ext,
							dmom_v_adv_w_ext,
							dmom_w_adv_p_ext,
							dmom_w_adv_u_ext,
							dmom_w_adv_v_ext,
							dmom_w_adv_w_ext,
							dflux_mass_u_ext,
							dflux_mass_v_ext,
							dflux_mass_w_ext,
							dflux_mom_u_adv_p_ext,
							dflux_mom_u_adv_u_ext,
							dflux_mom_u_adv_v_ext,
							dflux_mom_u_adv_w_ext,
							dflux_mom_v_adv_p_ext,
							dflux_mom_v_adv_u_ext,
							dflux_mom_v_adv_v_ext,
							dflux_mom_v_adv_w_ext,
							dflux_mom_w_adv_p_ext,
							dflux_mom_w_adv_u_ext,
							dflux_mom_w_adv_v_ext,
							dflux_mom_w_adv_w_ext);
	      //
	      //calculate the flux jacobian
	      //
	      ck.calculateGScale(G,normal,h_penalty);
	      h_penalty = 10.0/h_penalty;
	      //cek debug, do it the old way
	      h_penalty = 100.0/elementDiameter[eN];
	      for (int j=0;j<nDOF_trial_element;j++)
		{
		  register int j_nSpace = j*nSpace,ebN_local_kb_j=ebN_local_kb*nDOF_trial_element+j;
		  //cek debug
		  //ebqe_penalty_ext[ebNE_kb] = 10.0;
		  //
		  //cek todo add full stress on boundaries

		  fluxJacobian_p_p[j]=0.0;
		  fluxJacobian_p_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_p_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_p_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

		  fluxJacobian_u_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
							   ebqe_phi_ext[ebNE_kb],
							   sdInfo_u_u_rowptr,
							   sdInfo_u_u_colind,
							   isDOFBoundary_u[ebNE_kb],
							   normal,
							   mom_u_diff_ten_ext,
							   vel_trial_trace_ref[ebN_local_kb_j],
							   &vel_grad_trial_trace[j_nSpace],
							   h_penalty);//ebqe_penalty_ext[ebNE_kb]);
		  fluxJacobian_u_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_u_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

		  fluxJacobian_v_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_v_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_v_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
							   ebqe_phi_ext[ebNE_kb],
							   sdInfo_v_v_rowptr,
							   sdInfo_v_v_colind,
							   isDOFBoundary_v[ebNE_kb],
							   normal,
							   mom_v_diff_ten_ext,
							   vel_trial_trace_ref[ebN_local_kb_j],
							   &vel_grad_trial_trace[j_nSpace],
							   h_penalty);//ebqe_penalty_ext[ebNE_kb]);
		  fluxJacobian_v_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

		  fluxJacobian_w_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_w_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_w_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_w_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
							   ebqe_phi_ext[ebNE_kb],
							   sdInfo_w_w_rowptr,
							   sdInfo_w_w_colind,
							   isDOFBoundary_w[ebNE_kb],
							   normal,
							   mom_w_diff_ten_ext,
							   vel_trial_trace_ref[ebN_local_kb_j],
							   &vel_grad_trial_trace[j_nSpace],
							   h_penalty);//ebqe_penalty_ext[ebNE_kb]);

		  fluxJacobian_p_p[j]=0.0;
		  fluxJacobian_p_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_p_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_p_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mass_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

		  fluxJacobian_u_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_u_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
							   ebqe_phi_ext[ebNE_kb],
							   sdInfo_u_u_rowptr,
							   sdInfo_u_u_colind,
							   isDOFBoundary_u[ebNE_kb],
							   normal,
							   mom_u_diff_ten_ext,
							   vel_trial_trace_ref[ebN_local_kb_j],
							   &vel_grad_trial_trace[j_nSpace],
							   h_penalty);//ebqe_penalty_ext[ebNE_kb]);
		  fluxJacobian_u_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_u_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_u_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

		  fluxJacobian_v_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_v_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_v_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
							   ebqe_phi_ext[ebNE_kb],
							   sdInfo_v_v_rowptr,
							   sdInfo_v_v_colind,
							   isDOFBoundary_v[ebNE_kb],
							   normal,
							   mom_v_diff_ten_ext,
							   vel_trial_trace_ref[ebN_local_kb_j],
							   &vel_grad_trial_trace[j_nSpace],
							   h_penalty);//ebqe_penalty_ext[ebNE_kb]);
		  fluxJacobian_v_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_v_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]);

		  fluxJacobian_w_p[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_p_ext,p_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_w_u[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_u_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_w_v[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_v_ext,vel_trial_trace_ref[ebN_local_kb_j]);
		  fluxJacobian_w_w[j]=ck.ExteriorNumericalAdvectiveFluxJacobian(dflux_mom_w_adv_w_ext,vel_trial_trace_ref[ebN_local_kb_j]) +
		    ExteriorNumericalDiffusiveFluxJacobian(eps_rho,
							   ebqe_phi_ext[ebNE_kb],
							   sdInfo_w_w_rowptr,
							   sdInfo_w_w_colind,
							   isDOFBoundary_w[ebNE_kb],
							   normal,
							   mom_w_diff_ten_ext,
							   vel_trial_trace_ref[ebN_local_kb_j],
							   &vel_grad_trial_trace[j_nSpace],
							   h_penalty);//ebqe_penalty_ext[ebNE_kb]);
		  // //cek debug
		  // fluxJacobian_p_p[j]=0.0;
		  // fluxJacobian_p_u[j]=0.0;
		  // fluxJacobian_p_v[j]=0.0;
		  // fluxJacobian_p_w[j]=0.0;

		  // fluxJacobian_u_p[j]=0.0;
		  // fluxJacobian_u_u[j]=0.0;
		  // fluxJacobian_u_v[j]=0.0;
		  // fluxJacobian_u_w[j]=0.0;

		  // fluxJacobian_v_p[j]=0.0;
		  // fluxJacobian_v_u[j]=0.0;
		  // fluxJacobian_v_v[j]=0.0;
		  // fluxJacobian_v_w[j]=0.0;

		  // fluxJacobian_w_p[j]=0.0;
		  // fluxJacobian_w_u[j]=0.0;
		  // fluxJacobian_w_v[j]=0.0;
		  // fluxJacobian_w_w[j]=0.0;
		  // //cek debug
		}//j
	      //
	      //update the global Jacobian from the flux Jacobian
	      //
	      for (int i=0;i<nDOF_test_element;i++)
		{
		  register int eN_i = eN*nDOF_test_element+i;
		  for (int j=0;j<nDOF_trial_element;j++)
		    {
		      register int ebN_i_j = ebN*4*nDOF_test_X_trial_element + i*nDOF_trial_element + j;
		  
		      globalJacobian[csrRowIndeces_p_p[eN_i] + csrColumnOffsets_eb_p_p[ebN_i_j]] += fluxJacobian_p_p[j]*p_test_dS[i];
		      globalJacobian[csrRowIndeces_p_u[eN_i] + csrColumnOffsets_eb_p_u[ebN_i_j]] += fluxJacobian_p_u[j]*p_test_dS[i];
		      globalJacobian[csrRowIndeces_p_v[eN_i] + csrColumnOffsets_eb_p_v[ebN_i_j]] += fluxJacobian_p_v[j]*p_test_dS[i];
		      globalJacobian[csrRowIndeces_p_w[eN_i] + csrColumnOffsets_eb_p_w[ebN_i_j]] += fluxJacobian_p_w[j]*p_test_dS[i];
		   
		      globalJacobian[csrRowIndeces_u_p[eN_i] + csrColumnOffsets_eb_u_p[ebN_i_j]] += fluxJacobian_u_p[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_u_u[eN_i] + csrColumnOffsets_eb_u_u[ebN_i_j]] += fluxJacobian_u_u[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_u_v[eN_i] + csrColumnOffsets_eb_u_v[ebN_i_j]] += fluxJacobian_u_v[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_u_w[eN_i] + csrColumnOffsets_eb_u_w[ebN_i_j]] += fluxJacobian_u_w[j]*vel_test_dS[i];
		   
		      globalJacobian[csrRowIndeces_v_p[eN_i] + csrColumnOffsets_eb_v_p[ebN_i_j]] += fluxJacobian_v_p[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_v_u[eN_i] + csrColumnOffsets_eb_v_u[ebN_i_j]] += fluxJacobian_v_u[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_v_v[eN_i] + csrColumnOffsets_eb_v_v[ebN_i_j]] += fluxJacobian_v_v[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_v_w[eN_i] + csrColumnOffsets_eb_v_w[ebN_i_j]] += fluxJacobian_v_w[j]*vel_test_dS[i];
		   
		      globalJacobian[csrRowIndeces_w_p[eN_i] + csrColumnOffsets_eb_w_p[ebN_i_j]] += fluxJacobian_w_p[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_w_u[eN_i] + csrColumnOffsets_eb_w_u[ebN_i_j]] += fluxJacobian_w_u[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_w_v[eN_i] + csrColumnOffsets_eb_w_v[ebN_i_j]] += fluxJacobian_w_v[j]*vel_test_dS[i];
		      globalJacobian[csrRowIndeces_w_w[eN_i] + csrColumnOffsets_eb_w_w[ebN_i_j]] += fluxJacobian_w_w[j]*vel_test_dS[i];
		    }//j
		}//i
	      // //debug
	      // std::cout<<"flux jacobian ebNE "<<ebNE<<" kb "<<kb<<std::endl;
	      // for (int i=0;i<nDOF_test_element;i++)
	      //   {
	      //     for (int j=0;j<nDOF_trial_element;j++)
	      // 	{
	      // 	  std::cout<< fluxJacobian_p_p[j]*p_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_p_u[j]*p_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_p_v[j]*p_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_p_w[j]*p_test_dS[i]<<std::endl;
		  
	      // 	  std::cout<< fluxJacobian_u_p[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_u_u[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_u_v[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_u_w[j]*vel_test_dS[i]<<std::endl;
		  
	      // 	  std::cout<< fluxJacobian_v_p[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_v_u[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_v_v[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_v_w[j]*vel_test_dS[i]<<std::endl;
		  
	      // 	  std::cout<< fluxJacobian_w_p[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_w_u[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_w_v[j]*vel_test_dS[i]<<std::endl;
	      // 	  std::cout<< fluxJacobian_w_w[j]*vel_test_dS[i]<<std::endl;
	      // 	}//j
	      //   }//i
	    }//kb
	}//ebNE
    }//computeJacobian

    void calculateVelocityAverage(int nExteriorElementBoundaries_global,
    				  int* exteriorElementBoundariesArray,
    				  int nInteriorElementBoundaries_global,
    				  int* interiorElementBoundariesArray,
    				  int* elementBoundaryElementsArray,
    				  int* elementBoundaryLocalElementBoundariesArray,
    				  double* mesh_dof,
    				  int* mesh_l2g,
    				  double* mesh_trial_trace_ref,
    				  double* mesh_grad_trial_trace_ref,
    				  double* normal_ref,
    				  double* boundaryJac_ref,
    				  int* vel_l2g,
    				  double* u_dof,
    				  double* v_dof,
    				  double* w_dof,
    				  double* vel_trial_trace_ref,
    				  double* ebqe_velocity,
    				  double* velocityAverage)
    {
      int permutations[nQuadraturePoints_elementBoundary];
      double xArray_left[nQuadraturePoints_elementBoundary*3],
    	xArray_right[nQuadraturePoints_elementBoundary*3];
      for (int ebNE = 0; ebNE < nExteriorElementBoundaries_global; ebNE++)
    	{
    	  register int ebN = exteriorElementBoundariesArray[ebNE];
    	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
    	    {
    	      register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace,
    		ebNE_kb_nSpace = ebNE*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
    	      velocityAverage[ebN_kb_nSpace+0]=ebqe_velocity[ebNE_kb_nSpace+0];
    	      velocityAverage[ebN_kb_nSpace+1]=ebqe_velocity[ebNE_kb_nSpace+1];
    	      velocityAverage[ebN_kb_nSpace+2]=ebqe_velocity[ebNE_kb_nSpace+2];
    	    }//ebNE
    	}
      for (int ebNI = 0; ebNI < nInteriorElementBoundaries_global; ebNI++)
    	{
    	  register int ebN = interiorElementBoundariesArray[ebNI],
    	    left_eN_global   = elementBoundaryElementsArray[ebN*2+0],
    	    left_ebN_element  = elementBoundaryLocalElementBoundariesArray[ebN*2+0],
    	    right_eN_global  = elementBoundaryElementsArray[ebN*2+1],
    	    right_ebN_element = elementBoundaryLocalElementBoundariesArray[ebN*2+1],
    	    left_eN_nDOF_trial_element = left_eN_global*nDOF_trial_element,
    	    right_eN_nDOF_trial_element = right_eN_global*nDOF_trial_element;
    	  double jac[nSpace*nSpace],
    	    jacDet,
    	    jacInv[nSpace*nSpace],
    	    boundaryJac[nSpace*(nSpace-1)],
    	    metricTensor[(nSpace-1)*(nSpace-1)],
    	    metricTensorDetSqrt,
    	    normal[3],
    	    x,y,z;
    	  //double G[nSpace*nSpace],G_dd_G,tr_G,h_phi,h_penalty;
	  
    	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
    	    {
    	      ck.calculateMapping_elementBoundary(left_eN_global,
    						  left_ebN_element,
    						  kb,
    						  left_ebN_element*kb,
    						  mesh_dof,
    						  mesh_l2g,
    						  mesh_trial_trace_ref,
    						  mesh_grad_trial_trace_ref,
    						  boundaryJac_ref,
    						  jac,
    						  jacDet,
    						  jacInv,
    						  boundaryJac,
    						  metricTensor,
    						  metricTensorDetSqrt,
    						  normal_ref,
    						  normal,
    						  x,y,z);
    	      xArray_left[kb*3+0] = x;
    	      xArray_left[kb*3+1] = y;
    	      xArray_left[kb*3+2] = z;
    	      ck.calculateMapping_elementBoundary(right_eN_global,
    						  right_ebN_element,
    						  kb,
    						  right_ebN_element*kb,
    						  mesh_dof,
    						  mesh_l2g,
    						  mesh_trial_trace_ref,
    						  mesh_grad_trial_trace_ref,
    						  boundaryJac_ref,
    						  jac,
    						  jacDet,
    						  jacInv,
    						  boundaryJac,
    						  metricTensor,
    						  metricTensorDetSqrt,
    						  normal_ref,
    						  normal,
    						  x,y,z);
    	      xArray_right[kb*3+0] = x;
    	      xArray_right[kb*3+1] = y;
    	      xArray_right[kb*3+2] = z;
    	    }
    	  for  (int kb_left=0;kb_left<nQuadraturePoints_elementBoundary;kb_left++)
    	    {
    	      double errorNormMin = 1.0;
    	      for  (int kb_right=0;kb_right<nQuadraturePoints_elementBoundary;kb_right++)
    		{
    		  double errorNorm=0.0;
    		  for (int I=0;I<nSpace;I++)
    		    {
    		      errorNorm += fabs(xArray_left[kb_left*3+I]
    					-
    					xArray_right[kb_right*3+I]);
    		    }
    		  if (errorNorm < errorNormMin)
    		    {
    		      permutations[kb_right] = kb_left;
    		      errorNormMin = errorNorm;
    		    }
    		}
    	    }
    	  for  (int kb=0;kb<nQuadraturePoints_elementBoundary;kb++)
    	    {
    	      register int ebN_kb_nSpace = ebN*nQuadraturePoints_elementBoundary*nSpace+kb*nSpace;
    	      register double u_left=0.0,
    		v_left=0.0,
    		w_left=0.0,
    		u_right=0.0,
    		v_right=0.0,
    		w_right=0.0;
    	      register int left_kb = kb,
    		right_kb = permutations[kb],
    		left_ebN_element_kb_nDOF_test_element=left_ebN_element*left_kb*nDOF_test_element,
    		right_ebN_element_kb_nDOF_test_element=right_ebN_element*right_kb*nDOF_test_element;
    	      //
    	      //calculate the velocity solution at quadrature points on left and right
    	      //
    	      ck.valFromDOF(u_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],u_left);
    	      ck.valFromDOF(v_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],v_left);
    	      ck.valFromDOF(w_dof,&vel_l2g[left_eN_nDOF_trial_element],&vel_trial_trace_ref[left_ebN_element_kb_nDOF_test_element],w_left);
    	      //
    	      ck.valFromDOF(u_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],u_right);
    	      ck.valFromDOF(v_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],v_right);
    	      ck.valFromDOF(w_dof,&vel_l2g[right_eN_nDOF_trial_element],&vel_trial_trace_ref[right_ebN_element_kb_nDOF_test_element],w_right);
    	      //
    	      velocityAverage[ebN_kb_nSpace+0]=0.5*(u_left + u_right);
    	      velocityAverage[ebN_kb_nSpace+1]=0.5*(v_left + v_right);
    	      velocityAverage[ebN_kb_nSpace+2]=0.5*(w_left + w_right);
    	    }//ebNI
    	}
    }
  };//RANS3P
  
  inline RANS3P_base* newRANS3P(int nSpaceIn,
				int nQuadraturePoints_elementIn,
				int nDOF_mesh_trial_elementIn,
				int nDOF_trial_elementIn,
				int nDOF_test_elementIn,
				int nQuadraturePoints_elementBoundaryIn,
				int CompKernelFlag)
  {
    return proteus::chooseAndAllocateDiscretization<RANS3P_base,RANS3P,CompKernel>(nSpaceIn,
										   nQuadraturePoints_elementIn,
										   nDOF_mesh_trial_elementIn,
										   nDOF_trial_elementIn,
										   nDOF_test_elementIn,
										   nQuadraturePoints_elementBoundaryIn,
										   CompKernelFlag);
  }
}//proteus

#endif
