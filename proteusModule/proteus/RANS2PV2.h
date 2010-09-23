#ifndef RANS2PV2_H
#define RANS2PV2_H
#include <cmath>
#include <iostream>
#include "CompKernel.h"

#define nSpace 3
#define nQuadraturePoints_element 5
#define nDOF_trial_element 4
#define nDOF_mesh_trial_element 4
#define nDOF_test_element 4
#define nDOF_test_X_trial_element 16
#define nQuadraturePoints_elementBoundary 4
#define nElementBoundaries_element 4

namespace RANS2PV2
{
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
			      const double& phi,
			      const double n[nSpace],
			      const double& kappa,
			      const double& p,
			      const double grad_p[nSpace],
			      const double& u,
			      const double& v,
			      const double& w,
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
			      double& mom_v_source,
			      double& mom_w_source,
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
  
    //u momentum accumulation
    mom_u_acc=rho*u;
    dmom_u_acc_u=rho;
  
    //v momentum accumulation
    mom_v_acc=rho*v;
    dmom_v_acc_v=rho;
  
    //w momentum accumulation
    mom_w_acc=rho*w;
    dmom_w_acc_w=rho;
  
  
    //mass advective flux
    mass_adv[0]=u;
    mass_adv[1]=v;
    mass_adv[2]=w;
  
    dmass_adv_u[0]=1.0;
    dmass_adv_u[1]=0.0;
    dmass_adv_u[2]=0.0;

    dmass_adv_v[0]=0.0;
    dmass_adv_v[1]=1.0;
    dmass_adv_v[2]=0.0;

    dmass_adv_w[0]=0.0;
    dmass_adv_w[1]=0.0;
    dmass_adv_w[2]=1.0;

    //u momentum advective flux
    mom_u_adv[0]=rho*u*u;
    mom_u_adv[1]=rho*u*v;
    mom_u_adv[2]=rho*u*w;
  
    dmom_u_adv_u[0]=rho*2.0*u;
    dmom_u_adv_u[1]=rho*v;
    dmom_u_adv_u[2]=rho*w;
  
    dmom_u_adv_v[0]=0.0;
    dmom_u_adv_v[1]=rho*u;
    dmom_u_adv_v[2]=0.0;
  
    dmom_u_adv_w[0]=0.0;
    dmom_u_adv_w[1]=0.0;
    dmom_u_adv_w[2]=rho*u;
  
    //v momentum advective_flux
    mom_v_adv[0]=rho*v*u;
    mom_v_adv[1]=rho*v*v;
    mom_v_adv[2]=rho*v*w;
  
    dmom_v_adv_u[0]=rho*v;
    dmom_v_adv_u[1]=0.0;
    dmom_v_adv_u[2]=0.0;
  
    dmom_v_adv_w[0]=0.0;
    dmom_v_adv_w[1]=0.0;
    dmom_v_adv_w[2]=rho*v;
  
    dmom_v_adv_v[0]=rho*u;
    dmom_v_adv_v[1]=rho*2.0*v;
    dmom_v_adv_v[2]=rho*w;
  
    //w momentum advective_flux
    mom_w_adv[0]=rho*w*u;
    mom_w_adv[1]=rho*w*v;
    mom_w_adv[2]=rho*w*w;
  
    dmom_w_adv_u[0]=rho*w;
    dmom_w_adv_u[1]=0.0;
    dmom_w_adv_u[2]=0.0;
  
    dmom_w_adv_v[0]=0.0;
    dmom_w_adv_v[1]=rho*w;
    dmom_w_adv_v[2]=0.0;
  
    dmom_w_adv_w[0]=rho*u;
    dmom_w_adv_w[1]=rho*v;
    dmom_w_adv_w[2]=rho*2.0*w;
  
    //u momentum diffusion tensor
    mom_u_diff_ten[0] = 2.0*mu;
    mom_u_diff_ten[1] = mu;
    mom_u_diff_ten[2] = mu;
  
    mom_uv_diff_ten[0]=mu;
  
    mom_uw_diff_ten[0]=mu;
  
    //v momentum diffusion tensor
    mom_v_diff_ten[0] = mu;
    mom_v_diff_ten[1] = 2.0*mu;
    mom_v_diff_ten[2] = mu;
  
    mom_vu_diff_ten[0]=mu;
  
    mom_vw_diff_ten[0]=mu;
  
    //w momentum diffusion tensor
    mom_w_diff_ten[0] = mu;
    mom_w_diff_ten[1] = mu;
    mom_w_diff_ten[2] = 2.0*mu;
  
    mom_wu_diff_ten[0]=mu;
  
    mom_wv_diff_ten[0]=mu;
  
    //momentum sources
    norm_n = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    mom_u_source = -rho*g[0] - d_mu*sigma*kappa*n[0];
    mom_v_source = -rho*g[1] - d_mu*sigma*kappa*n[1];
    mom_w_source = -rho*g[2] - d_mu*sigma*kappa*n[2];
   
    //u momentum Hamiltonian (pressure)
    mom_u_ham = grad_p[0];
    dmom_u_ham_grad_p[0]=1.0;
    dmom_u_ham_grad_p[1]=0.0;
    dmom_u_ham_grad_p[2]=0.0;
  
    //v momentum Hamiltonian (pressure)
    mom_v_ham = grad_p[1];
    dmom_v_ham_grad_p[0]=0.0;
    dmom_v_ham_grad_p[1]=1.0;
    dmom_v_ham_grad_p[2]=0.0;
  
    //w momentum Hamiltonian (pressure)
    mom_w_ham = grad_p[2];
    dmom_w_ham_grad_p[0]=0.0;
    dmom_w_ham_grad_p[1]=0.0;
    dmom_w_ham_grad_p[2]=1.0;
  }
  

  inline
    void calculateSubgridError_tau(const double&  Ct_sge,
				   const double&  Cd_sge,
				   const double* G,
				   const double& G_dd_G,
				   const double& tr_G,
				   const double& rho,
				   const double& Dt,
				   const double v[nSpace],
				   const double& mu,
				   double& tau_p,
				   double& tau_v,
				   double& cfl)
  {
    const double rho2=rho*rho,Dt2=Dt*Dt,mu2=mu*mu;
    register double v_d_Gv=0.0;
    for(int I=0;I<nSpace;I++)
      for (int J=0;J<nSpace;J++)
	v_d_Gv += v[I]*G[I*nSpace+J]*v[J];
    cfl = 2.0*sqrt(v_d_Gv);
    //cek 1.0/sqrt(rho2*Dt2 + 4*v_d_Gv + 144*mu2*G_dd_G); ?
    /* tau_v = 1.0/sqrt(rho2*Dt2 + 4*v_d_Gv + 144*mu2*G_dd_G); */
    /* tau_p = 1.0/(tr_G*tau_v);  */
    //cek "correct" tau
    tau_v = 1.0/sqrt(Ct_sge*rho2*Dt2 + rho2*v_d_Gv + Cd_sge*mu2*G_dd_G);
    tau_p = 1.0/(tr_G*tau_v);
    //debug
    /* double tau_v_old = tau_v,tau_p_old = tau_p; */
    /* double nrm_v=0.0,h=1.0/20.0; */
    /* double oneByAbsdt =  fabs(Dt); */
    /* for(int I=0;I<nSpace;I++) */
    /* 	nrm_v += v[I]*v[I]; */
    /* nrm_v = sqrt(nrm_v); */
    /* cfl = nrm_v/h; */
    /* tau_v = 1.0/(4.0*mu/(h*h) + 2.0*rho*nrm_v/h + oneByAbsdt); */
    /* tau_p = 4.0*mu + 2.0*rho*nrm_v*h+ oneByAbsdt*h*h; */
    /* std::cout<<"nrm_v "<<nrm_v<<" tau_v "<<tau_v<<"\t"<<tau_v_old<<" tau_p "<<tau_p<<'\t'<<tau_p_old<<std::endl; */
  }

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
	flux_umom+= n[0]*(bc_p-p);
	flux_vmom+= n[1]*(bc_p-p);
	flux_wmom+= n[2]*(bc_p-p);
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
	dflux_umom_dp= -n[0];
	dflux_vmom_dp= -n[1];
	dflux_wmom_dp= -n[2];
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
	penaltyFlux = penalty*(u-bc_u);
	flux += penaltyFlux;
	//contact line slip
	flux*=(smoothedDirac(eps,0) - smoothedDirac(eps,phi))/smoothedDirac(eps,0);
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
    double dvel_I,tmp=0.0;
    if(isDOFBoundary >= 1)
      {
	for(int I=0;I<nSpace;I++)
	  {
	    dvel_I=0.0;
	    for(int m=rowptr[I];m<rowptr[I+1];m++)
	      {
		dvel_I -= a[m]*grad_v[colind[m]];
	      }
	    tmp += dvel_I*n[I];
	  }
	tmp +=penalty*v;
	//contact line slip
	tmp*=(smoothedDirac(eps,0) - smoothedDirac(eps,phi))/smoothedDirac(eps,0);
      }
    return tmp;
  }

}//RANS2PV2
extern "C"
{
  void calculateResidual_RANS2PV2(//testing mesh replacement
				  double* mesh_trial_ref,
				  double* mesh_grad_trial_ref,
				  double* mesh_dof,
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
				  //end testing meshreplacement
				  int nElements_global,
				  double alpha_bdf,
				  double eps_rho,
				  double eps_mu,
				  double sigma,
				  double rho_0,
				  double nu_0,
				  double rho_1,
				  double nu_1,
				  double Ct_sge,
				  double Cd_sge,
				  double C_dc,
				  int* p_l2g, int* vel_l2g,
				  double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				  double* g,
				  double* phi,
				  double* n,
				  double* kappa,
				  double* q_mom_u_acc,
				  double* q_mom_v_acc,
				  double* q_mom_w_acc,
				  double* q_mass_adv,
				  double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				  double* q_velocity_last,
				  double* q_cfl,
				  double* q_numDiff_u, double* q_numDiff_v, double* q_numDiff_w,
				  double* q_numDiff_u_last, double* q_numDiff_v_last, double* q_numDiff_w_last,
				  double* q_elementResidual_p, double* q_elementResidual_u, double* q_elementResidual_v, double* q_elementResidual_w,
				  int* sdInfo_u_u_rowptr,int* sdInfo_u_u_colind,			      
				  int* sdInfo_u_v_rowptr,int* sdInfo_u_v_colind,
				  int* sdInfo_u_w_rowptr,int* sdInfo_u_w_colind,
				  int* sdInfo_v_v_rowptr,int* sdInfo_v_v_colind,
				  int* sdInfo_v_u_rowptr,int* sdInfo_v_u_colind,
				  int* sdInfo_v_w_rowptr,int* sdInfo_v_w_colind,
				  int* sdInfo_w_w_rowptr,int* sdInfo_w_w_colind,
				  int* sdInfo_w_u_rowptr,int* sdInfo_w_u_colind,
				  int* sdInfo_w_v_rowptr,int* sdInfo_w_v_colind,
				  int offset_p, int offset_u, int offset_v, int offset_w, int stride_p, int stride_u, int stride_v, int stride_w, double* globalResidual,
				  int nExteriorElementBoundaries_global,
				  int* exteriorElementBoundariesArray,
				  int* elementBoundaryElementsArray,
				  int* elementBoundaryLocalElementBoundariesArray,
				  double* ebqe_phi_ext,
				  double* ebqe_n_ext,
				  double* ebqe_kappa_ext,
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
				  double* ebqe_velocity_ext,
				  double* flux);

  void calculateJacobian_RANS2PV2(//testing mesh replacement
				  double* mesh_trial_ref,
				  double* mesh_grad_trial_ref,
				  double* mesh_dof,
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
				  //end testing meshreplacement
				  int nElements_global,
				  double alpha_bdf,
				  double eps_rho,
				  double eps_mu,
				  double sigma,
				  double rho_0,
				  double nu_0,
				  double rho_1,
				  double nu_1,
				  double Ct_sge,
				  double Cd_sge,
				  double C_dc,
				  int* p_l2g, int* vel_l2g,
				  double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				  double* g,
				  double* phi,
				  double* n,
				  double* kappa,
				  double* q_mom_u_acc_beta_bdf, double* q_mom_v_acc_beta_bdf, double* q_mom_w_acc_beta_bdf,
				  double* q_velocity_last,
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
				  double* ebqe_n_ext,
				  double* ebqe_kappa_ext,
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
				  int* csrColumnOffsets_eb_w_w);
  void calculateVelocityAverage_RANS2PV2(int* permutations,
					 int nExteriorElementBoundaries_global,
					 int* exteriorElementBoundariesArray,
					 int nInteriorElementBoundaries_global,
					 int* interiorElementBoundariesArray,
					 int* elementBoundaryElementsArray,
					 int* elementBoundaryLocalElementBoundariesArray,
					 int* vel_l2g, 
					 double* u_dof, double* v_dof, double* w_dof,
					 double* vel_trial,
					 double* ebqe_velocity,
					 double* velocityAverage);
}//extern "C"
#endif
