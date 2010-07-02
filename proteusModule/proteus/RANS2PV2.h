#ifndef RANS2PV2_H
#define RANS2PV2_H
#include <cmath>
#include <iostream>

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
inline
double smoothedHeaviside(double eps, double phi)
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

inline
double smoothedHeaviside_integral(double eps, double phi)
{
  double HI;
  if (phi > eps)
    {
      HI= phi - eps + 	0.5*(eps + 0.5*eps*eps/eps - eps*cos(M_PI*eps/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  else if (phi < -eps)
    {
      HI=0.0;
    }
  else
    {
      HI = 0.5*(phi + 0.5*phi*phi/eps - eps*cos(M_PI*phi/eps)/(M_PI*M_PI)) - 0.5*((-eps) + 0.5*(-eps)*(-eps)/eps - eps*cos(M_PI*(-eps)/eps)/(M_PI*M_PI));
    }
  return HI;
}

inline
double smoothedDirac(double eps, double phi)
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
double valFromDOF_c(const double& dof,const double& v)
{
  return dof*v;
}

inline
void valFromDOF(int eN, int k, int j, double* u,double* dof,int* l2g,double* v)
{
  u[eN*nQuadraturePoints_element+
    k]
    += valFromDOF_c(dof[l2g[eN*nDOF_trial_element+
			    j]],
		    v[eN*nQuadraturePoints_element*nDOF_trial_element+
		      k*nDOF_trial_element+
		      j]);
}


inline
void valFromDOF_ext(int ebNE, int eN, int k, int j, double* u, double* dof, int* l2g, double* v)
{
  u[ebNE*nQuadraturePoints_elementBoundary+
    k] 
    += valFromDOF_c(dof[l2g[eN*nDOF_trial_element+j]],
		    v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+
		      k*nDOF_trial_element+
		      j]);
}


inline
double gradFromDOF_c(const double& dof,const double& grad_v_I)
{
  return dof*grad_v_I;
}


inline
void gradFromDOF(int eN, int k, int j, double* grad_u,double* dof,int* l2g,double* grad_v)
{
  for (int I=0;I<nSpace;I++)
    grad_u[eN*nQuadraturePoints_element*nSpace+
	   k*nSpace+
	   I] += gradFromDOF_c(dof[l2g[eN*nDOF_trial_element+
				       j]],
			       grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace+
				      k*nDOF_trial_element*nSpace+
				      j*nSpace+
				      I]);
}

inline
void gradFromDOF_ext(int ebNE, int eN, int k, int j, double* grad_u, double* dof, int* l2g, double* grad_v)
{
  for (int I=0;I<nSpace;I++)
    grad_u[ebNE*nQuadraturePoints_elementBoundary*nSpace+
	   k*nSpace+I]
      += gradFromDOF_c(dof[l2g[eN*nDOF_trial_element+j]],
		       grad_v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element*nSpace+
			      k*nDOF_trial_element*nSpace+
			      j*nSpace+I]);
}

inline
void evaluateCoefficients_c(const double eps_rho,
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
  mom_u_acc=u;
  dmom_u_acc_u=1.0;
  
  //v momentum accumulation
  mom_v_acc=v;
  dmom_v_acc_v=1.0;
  
  //w momentum accumulation
  mom_w_acc=w;
  dmom_w_acc_w=1.0;
  
  
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
  mom_u_adv[0]=u*u;
  mom_u_adv[1]=u*v;
  mom_u_adv[2]=u*w;
  
  dmom_u_adv_u[0]=2.0*u;
  dmom_u_adv_u[1]=v;
  dmom_u_adv_u[2]=w;
  
  dmom_u_adv_v[0]=0.0;
  dmom_u_adv_v[1]=u;
  dmom_u_adv_v[2]=0.0;
  
  dmom_u_adv_w[0]=0.0;
  dmom_u_adv_w[1]=0.0;
  dmom_u_adv_w[2]=u;
  
  //v momentum advective_flux
  mom_v_adv[0]=v*u;
  mom_v_adv[1]=v*v;
  mom_v_adv[2]=v*w;
  
  dmom_v_adv_u[0]=v;
  dmom_v_adv_u[1]=0.0;
  dmom_v_adv_u[2]=0.0;
  
  dmom_v_adv_w[0]=0.0;
  dmom_v_adv_w[1]=0.0;
  dmom_v_adv_w[2]=v;
  
  dmom_v_adv_v[0]=u;
  dmom_v_adv_v[1]=2.0*v;
  dmom_v_adv_v[2]=w;
  
  //w momentum advective_flux
  mom_w_adv[0]=w*u;
  mom_w_adv[1]=w*v;
  mom_w_adv[2]=w*w;
  
  dmom_w_adv_u[0]=w;
  dmom_w_adv_u[1]=0.0;
  dmom_w_adv_u[2]=0.0;
  
  dmom_w_adv_v[0]=0.0;
  dmom_w_adv_v[1]=w;
  dmom_w_adv_v[2]=0.0;
  
  dmom_w_adv_w[0]=u;
  dmom_w_adv_w[1]=v;
  dmom_w_adv_w[2]=2.0*w;
  
  //u momentum diffusion tensor
  mom_u_diff_ten[0] = 2.0*nu;
  mom_u_diff_ten[1] = nu;
  mom_u_diff_ten[2] = nu;
  
  mom_uv_diff_ten[0]=nu;
  
  mom_uw_diff_ten[0]=nu;
  
  //v momentum diffusion tensor
  mom_v_diff_ten[0] = nu;
  mom_v_diff_ten[1] = 2.0*nu;
  mom_v_diff_ten[2] = nu;
  
  mom_vu_diff_ten[0]=nu;
  
  mom_vw_diff_ten[0]=nu;
  
  //w momentum diffusion tensor
  mom_w_diff_ten[0] = nu;
  mom_w_diff_ten[1] = nu;
  mom_w_diff_ten[2] = 2.0*nu;
  
  mom_wu_diff_ten[0]=nu;
  
  mom_wv_diff_ten[0]=nu;
  
  //momentum sources
  norm_n = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  mom_u_source = -g[0] - d_mu*sigma*kappa*n[0]/(rho*(norm_n+1.0e-8));
  mom_v_source = -g[1] - d_mu*sigma*kappa*n[1]/(rho*(norm_n+1.0e-8));
  mom_w_source = -g[2] - d_mu*sigma*kappa*n[2]/(rho*(norm_n+1.0e-8));
   
  //u momentum Hamiltonian (pressure)
  mom_u_ham = grad_p[0]/rho;
  dmom_u_ham_grad_p[0]=1.0/rho;
  dmom_u_ham_grad_p[1]=0.0;
  dmom_u_ham_grad_p[2]=0.0;
  
  //v momentum Hamiltonian (pressure)
  mom_v_ham = grad_p[1]/rho;
  dmom_v_ham_grad_p[0]=0.0;
  dmom_v_ham_grad_p[1]=1.0/rho;
  dmom_v_ham_grad_p[2]=0.0;
  
  //w momentum Hamiltonian (pressure)
  mom_w_ham = grad_p[2]/rho;
  dmom_w_ham_grad_p[0]=0.0;
  dmom_w_ham_grad_p[1]=0.0;
  dmom_w_ham_grad_p[2]=1.0/rho;
}
  
inline
void evaluateCoefficients(int K,
                          const double eps_rho,
                          const double eps_mu,
                          const double sigma,
                          const double rho_0,
                          const double nu_0,
                          const double rho_1,
                          const double nu_1,
                          const double* g,
                          const double* phi,
                          const double* n,
                          const double* kappa,
                          const double *p,
                          const double *grad_p,
                          const double *u,
                          const double *v,
                          const double *w,
                          double *mom_u_acc,
                          double *dmom_u_acc_u,
                          double *mom_v_acc,
                          double *dmom_v_acc_v,
                          double *mom_w_acc,
                          double *dmom_w_acc_w,
                          double *mass_adv,
                          double *dmass_adv_u,
                          double *dmass_adv_v,
                          double *dmass_adv_w,
                          double *mom_u_adv,
                          double *dmom_u_adv_u,
                          double *dmom_u_adv_v,
                          double *dmom_u_adv_w,
                          double *mom_v_adv,
                          double *dmom_v_adv_u,
                          double *dmom_v_adv_v,
                          double *dmom_v_adv_w,
                          double *mom_w_adv,
                          double *dmom_w_adv_u,
                          double *dmom_w_adv_v,
                          double *dmom_w_adv_w,
                          double *mom_u_diff_ten,
                          double *mom_v_diff_ten,
                          double *mom_w_diff_ten,
                          double *mom_uv_diff_ten,
                          double *mom_uw_diff_ten,
                          double *mom_vu_diff_ten,
                          double *mom_vw_diff_ten,
                          double *mom_wu_diff_ten,
                          double *mom_wv_diff_ten,
                          double *mom_u_source,
                          double *mom_v_source,
                          double *mom_w_source,
                          double *mom_u_ham,
                          double *dmom_u_ham_grad_p,
                          double *mom_v_ham,
                          double *dmom_v_ham_grad_p,
                          double *mom_w_ham,
                          double *dmom_w_ham_grad_p)
{
  evaluateCoefficients_c(eps_rho,
                         eps_mu,
                         sigma,
                         rho_0,
                         nu_0,
                         rho_1,
                         nu_1,
                         g,
                         phi[K],
                         &n[K*nSpace],
                         kappa[K],
                         p[K],
                         &grad_p[K*nSpace],
                         u[K],
                         v[K],
                         w[K],
                         mom_u_acc[K],
                         dmom_u_acc_u[K],
                         mom_v_acc[K],
                         dmom_v_acc_v[K],
                         mom_w_acc[K],
                         dmom_w_acc_w[K],
                         &mass_adv[K*nSpace],
                         &dmass_adv_u[K*nSpace],
                         &dmass_adv_v[K*nSpace],
                         &dmass_adv_w[K*nSpace],
                         &mom_u_adv[K*nSpace],
                         &dmom_u_adv_u[K*nSpace],
                         &dmom_u_adv_v[K*nSpace],
                         &dmom_u_adv_w[K*nSpace],
                         &mom_v_adv[K*nSpace],
                         &dmom_v_adv_u[K*nSpace],
                         &dmom_v_adv_v[K*nSpace],
                         &dmom_v_adv_w[K*nSpace],
                         &mom_w_adv[K*nSpace],
                         &dmom_w_adv_u[K*nSpace],
                         &dmom_w_adv_v[K*nSpace],
                         &dmom_w_adv_w[K*nSpace],
                         &mom_u_diff_ten[K*nSpace],
                         &mom_v_diff_ten[K*nSpace],
                         &mom_w_diff_ten[K*nSpace],
                         &mom_uv_diff_ten[K],
                         &mom_uw_diff_ten[K],
                         &mom_vu_diff_ten[K],
                         &mom_vw_diff_ten[K],
                         &mom_wu_diff_ten[K],
                         &mom_wv_diff_ten[K],
                         mom_u_source[K],
                         mom_v_source[K],
                         mom_w_source[K],
                         mom_u_ham[K],
                         &dmom_u_ham_grad_p[K*nSpace],
                         mom_v_ham[K],
                         &dmom_v_ham_grad_p[K*nSpace],
                         mom_w_ham[K],
                         &dmom_w_ham_grad_p[K*nSpace]);
}
  
inline
void backwardEuler_c(const double& dt, const double& m_old, const double& m, const double& dm, double& mt, double& dmt)
{  
  mt =(m-m_old)/dt;
  dmt = dm/dt;
}

inline
void backwardEuler(int eN, int k, double dt, double* m_old, double* m, double* dm, double* mt, double* dmt)
{  
  backwardEuler_c(dt,
		  m_old[eN*nQuadraturePoints_element+
			k],
		  m[eN*nQuadraturePoints_element+
		    k],
		  dm[eN*nQuadraturePoints_element+
		     k],
		  mt[eN*nQuadraturePoints_element+
		     k], 
		  dmt[eN*nQuadraturePoints_element+
		      k]);
}
inline
void bdf_c(const double& alpha, const double& beta, const double& m, const double& dm, double& mt, double& dmt)
{  
  mt =alpha*m + beta;
  dmt = alpha*dm;
}
inline
void bdf(int eN, int k, double alpha, double* beta, double* m, double* dm, double* mt, double* dmt)
{  
  bdf_c(alpha,
	beta[eN*nQuadraturePoints_element+
	     k],
	m[eN*nQuadraturePoints_element+
	  k],
	dm[eN*nQuadraturePoints_element+
	   k],
	mt[eN*nQuadraturePoints_element+
	   k], 
	dmt[eN*nQuadraturePoints_element+
	    k]);
}

inline
double Mass_weak_c(const double& mt,
		   const double& w_dV)
{
  return mt*w_dV;
}

inline
void updateMass_weak(int eN,
                     int k,
                     int i,
                     double* mt,
                     double* w_dV,
                     double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + 
		i]  += Mass_weak_c(mt[eN*nQuadraturePoints_element + 
				      k],
				   w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
					k*nDOF_test_element + 
					i]);
}

inline
double MassJacobian_weak_c(const double& dmt,
			   const double& v,
			   const double& w_dV)
{
  return dmt*v*w_dV;
}

inline
void updateMassJacobian_weak(int eN,
			     int k,
			     int j,
			     int i,
			     double* dmt,
			     double* v,
			     double* w_dV,
			     double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element+
			 j] += MassJacobian_weak_c(dmt[eN*nQuadraturePoints_element + 
						       k],
						   v[eN*nQuadraturePoints_element*nDOF_trial_element + 
						     k*nDOF_trial_element + 
						     j],
						   w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
							k*nDOF_test_element + 
							i]);
}

inline
double Mass_strong_c(const double& mt)
{
  return mt; 
}

inline
void updateMass_strong(int eN,
		       int k,
		       double* mt,
		       double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] += Mass_strong_c(mt[eN*nQuadraturePoints_element+
					 k]);
}

inline
double MassJacobian_strong_c(const double& dmt,
			     const double& v)
{
  return dmt*v;
}

inline
void updateMassJacobian_strong(int eN,
                               int k,
                               int j,
                               double* dmt,
                               double* v,
                               double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += MassJacobian_strong_c(dmt[eN*nQuadraturePoints_element+
						   k],
					       v[eN*nQuadraturePoints_element*nDOF_trial_element+
						 k*nDOF_trial_element + 
						 j]);
}

inline
double Mass_adjoint_c(const double& dmt,
		      const double& w_dV)
{
  return dmt*w_dV;
}

inline
void updateMass_adjoint(int eN,
			int k,
			int i,
			double* dmt,
			double* w_dV,
			double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Mass_adjoint_c(dmt[eN*nQuadraturePoints_element + 
				      k],
				  w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
				       k*nDOF_test_element + 
				       i]);
}

inline
double Advection_weak_c(const double  f[nSpace],
			const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= f[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateAdvection_weak(int eN,
                          int k,
                          int i,
                          double* f,
                          double* grad_w_dV,
                          double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + 
		i] += Advection_weak_c(&f[eN*nQuadraturePoints_element*nSpace + 
					  k*nSpace],
				       &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
						  k*nDOF_test_element*nSpace + 
						  i*nSpace]);
}

inline
double AdvectionJacobian_weak_c(const double df[nSpace],
				const double& v,
				const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= df[I]*v*grad_w_dV[I];
  return tmp;
}

inline
void updateAdvectionJacobian_weak(int eN,
				  int k,
				  int j,
				  int i,
				  double* df,
				  double* v,
				  double* grad_w_dV,
				  double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element+
			 j] += AdvectionJacobian_weak_c(&df[eN*nQuadraturePoints_element*nSpace + 
							    k*nSpace],
							v[eN*nQuadraturePoints_element*nDOF_trial_element + 
							  k*nDOF_trial_element + 
							  j],
							&grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
								   k*nDOF_test_element*nSpace + 
								   i*nSpace]);
}

inline
double Advection_strong_c(const double df[nSpace],
			  const double grad_u[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += df[I]*grad_u[I];
  return tmp;
}

inline
void updateAdvection_strong(int eN,
			    int k,
			    double* df,
			    double* grad_u,
			    double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] += Advection_strong_c(&df[eN*nQuadraturePoints_element*nSpace + 
					       k*nSpace],
					   &grad_u[eN*nQuadraturePoints_element*nSpace + 
						   k*nSpace]);
}

inline
double AdvectionJacobian_strong_c(const double df[nSpace],
				  const double grad_v[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += df[I]*grad_v[I];
  return tmp;
}

inline
void updateAdvectionJacobian_strong(int eN,
                                    int k,
                                    int j,
                                    double* df,
                                    double* grad_v,
                                    double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += AdvectionJacobian_strong_c(&df[eN*nQuadraturePoints_element*nSpace + 
							k*nSpace],
						    &grad_v[eN*nQuadraturePoints_element*nSpace*nDOF_trial_element + 
							    k*nSpace*nDOF_trial_element +
							    j*nSpace]);
}

inline
double Advection_adjoint_c(const double df[nSpace],
			   const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= df[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateAdvection_adjoint(int eN,
			     int k,
			     int i,
			     double* df,
			     double* grad_w_dV,
			     double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Advection_adjoint_c(&df[eN*nQuadraturePoints_element*nSpace + 
					   k*nSpace],
				       &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
						  k*nDOF_test_element*nSpace + 
						  i*nSpace]);
}

inline
double Hamiltonian_weak_c(const double& H,
			  const double& w_dV)
{
  return H*w_dV;
}

inline
void updateHamiltonian_weak(int eN,
                            int k,
                            int i,
                            double* H,
                            double* w_dV,
                            double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] += Hamiltonian_weak_c(H[eN*nQuadraturePoints_element + 
								  k],
								w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
								     k*nDOF_test_element + 
								     i]);
}

inline
double HamiltonianJacobian_weak_c(const double dH[nSpace],
				  const double grad_v[nSpace],
				  const double& w_dV)
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += dH[I]*grad_v[I]*w_dV;
  return tmp;
}

inline
void updateHamiltonianJacobian_weak(int eN,
				    int k,
				    int j,
				    int i,
				    double* dH,
				    double* grad_v,
				    double* w_dV,
				    double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element +
			 j]+= HamiltonianJacobian_weak_c(&dH[eN*nQuadraturePoints_element*nSpace + 
							     k*nSpace],
							 &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace + 
								 k*nDOF_trial_element*nSpace +
								 j*nSpace],
							 w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
							      k*nDOF_test_element +
							      i]);
}

inline
double Hamiltonian_strong_c(const double dH[nSpace],
			    const double grad_u[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += dH[I]*grad_u[I];
  return tmp;
}

inline
void updateHamiltonian_strong(int eN,
			      int k,
			      double* dH,
			      double* grad_u,
			      double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+k] += Hamiltonian_strong_c(&dH[eN*nQuadraturePoints_element*nSpace+
									      k*nSpace],
									  &grad_u[eN*nQuadraturePoints_element*nSpace+
										  k*nSpace]);
}

inline
double HamiltonianJacobian_strong_c(const double dH[nSpace],
				    const double grad_v[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp += dH[I]*grad_v[I];
  return tmp;
}

inline
void updateHamiltonianJacobian_strong(int eN,
                                      int k,
                                      int j,
                                      double* dH,
                                      double* grad_v,
                                      double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += HamiltonianJacobian_strong_c(&dH[eN*nQuadraturePoints_element*nSpace+
							  k*nSpace],
						      &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace+
							      k*nDOF_trial_element*nSpace + 
							      j*nSpace]);
}

inline
double Hamiltonian_adjoint_c(const double dH[nSpace],
			     const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    tmp -= dH[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateHamiltonian_adjoint(int eN,
			       int k,
			       int i,
			       double* dH,
			       double* grad_w_dV,
			       double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Hamiltonian_adjoint_c(&dH[eN*nQuadraturePoints_element*nSpace + 
					     k*nSpace],
					 &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
						    k*nDOF_test_element*nSpace + 
						    i*nSpace]);
}

inline
double Diffusion_weak_c(int* rowptr,
			int* colind,
			double* a,
			const double grad_phi[nSpace],
			const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for(int I=0;I<nSpace;I++)
    for (int m=rowptr[I];m<rowptr[I+1];m++)
      tmp += a[m]*grad_phi[colind[m]]*grad_w_dV[I];
  return tmp;
}

inline
void updateDiffusion_weak(int eN,
			  int k,
			  int i,
			  int* rowptr,
			  int* colind,
			  double* a,
			  double* grad_phi,
			  double* grad_w_dV,
			  double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] += Diffusion_weak_c(rowptr,
							      colind,
							      &a[eN*nQuadraturePoints_element*rowptr[nSpace]+
								 k*rowptr[nSpace]],
							      &grad_phi[eN*nQuadraturePoints_element*nSpace + 
									k*nSpace],
							      &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
									 k*nDOF_test_element*nSpace + 
									 i*nSpace]);
}

inline
double DiffusionJacobian_weak_c(int* rowptr,
				int* colind,
				double* a,
				double* da,
				const double grad_phi[nSpace],
				const double grad_w_dV[nSpace],
				const double& dphi,
				const double& v,
				const double grad_v[nSpace])
{
  double daProduct=0.0,dphiProduct=0.0;
  for (int I=0;I<nSpace;I++)
    for (int m=rowptr[I];m<rowptr[I+1];m++)
      {
	daProduct += da[m]*grad_phi[colind[m]]*grad_w_dV[I];
	dphiProduct += a[m]*grad_v[colind[m]]*grad_w_dV[I];
      }
  return daProduct*v+dphiProduct*dphi;
}

inline
double SimpleDiffusionJacobian_weak_c(int* rowptr,
				      int* colind,
				      double* a,
				      const double grad_v[nSpace],
				      const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    for(int m=rowptr[I];m<rowptr[I+1];m++)
      tmp += a[m]*grad_v[colind[m]]*grad_w_dV[I];
  return tmp;
}

inline
void updateDiffusionJacobian_weak(int eN,
				  int k,
				  int j,
				  int i,
				  int* rowptr,
				  int* colind,
				  int* l2g,
				  double* a,
				  double* da,
				  double* grad_phi,
				  double* grad_w_dV,
				  double* dphi,
				  double* v,
				  double* grad_v,
				  double* jacobian_weak_residual)
{
  const int nnz=rowptr[nSpace];
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element + 
			 j] +=
    DiffusionJacobian_weak_c(rowptr,
			     colind,
			     &a[eN*nQuadraturePoints_element*nnz + 
				k*nnz],
			     &da[eN*nQuadraturePoints_element*nnz+
				 k*nnz],
			     &grad_phi[eN*nQuadraturePoints_element*nSpace +
				       k*nSpace],
			     &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace +
					k*nDOF_test_element*nSpace +
					i*nSpace],
			     dphi[l2g[eN*nDOF_trial_element + 
				      j]],
			     v[eN*nQuadraturePoints_element*nDOF_trial_element+
			       k*nDOF_trial_element+
			       j],
			     &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace +
				     k*nDOF_trial_element*nSpace +
				     j*nSpace]);
}

inline
double Reaction_weak_c(const double& r,
		       const double& w_dV)
{
  return r*w_dV;
}

inline
void updateReaction_weak(int eN,
			 int k,
			 int i,
			 double* r,
			 double* w_dV,
			 double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + 
		i]  += Reaction_weak_c(r[eN*nQuadraturePoints_element + 
					 k],
				       w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
					    k*nDOF_test_element + 
					    i]);
}

inline
double ReactionJacobian_weak_c(const double& dr,
			       const double& v,
			       const double& w_dV)
{
  return dr*v*w_dV;
}

inline
void updateReactionJacobian_weak(int eN,
				 int k,
				 int j,
				 int i,
				 double* dr,
				 double* v,
				 double* w_dV,
				 double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element+
			 j] += ReactionJacobian_weak_c(dr[eN*nQuadraturePoints_element + 
							  k],
						       v[eN*nQuadraturePoints_element*nDOF_trial_element + 
							 k*nDOF_trial_element + 
							 j],
						       w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
							    k*nDOF_test_element + 
							    i]);
}

inline
double Reaction_strong_c(const double& r)
{
  return r; 
}

inline
void updateReaction_strong(int eN,
			   int k,
			   double* r,
			   double* strong_residual)
{
  strong_residual[eN*nQuadraturePoints_element+
		  k] += Reaction_strong_c(r[eN*nQuadraturePoints_element+
					    k]);
}

inline
double ReactionJacobian_strong_c(const double& dr,
				     const double& v)
{
  return dr*v;
}

inline
void updateReactionJacobian_strong(int eN,
				   int k,
				   int j,
				   double* dr,
				   double* v,
				   double* dstrong_residual)
{
  dstrong_residual[eN*nQuadraturePoints_element*nDOF_trial_element+
		   k*nDOF_trial_element + 
		   j] += ReactionJacobian_strong_c(dr[eN*nQuadraturePoints_element+
						      k],
						   v[eN*nQuadraturePoints_element*nDOF_trial_element+
						     k*nDOF_trial_element + 
						     j]);
}

inline
double Reaction_adjoint_c(const double& dr,
			  const double& w_dV)
{
  return dr*w_dV;
}

inline
void updateReaction_adjoint(int eN,
			    int k,
			    int i,
			    double* dr,
			    double* w_dV,
			    double* Lstar_w_dV)
{
  Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
	     k*nDOF_test_element + 
	     i] += Reaction_adjoint_c(dr[eN*nQuadraturePoints_element + 
					 k],
				      w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
					   k*nDOF_test_element + 
					   i]);
}

inline
void calculateSubgridError_tau_c(const double&  hFactor,
				 const double& elementDiameter,
				 const double& dmt,
				 const double& dm,
				 const double f[nSpace],
				 const double& a,
				 double& tau0,
				 double& tau1,
				 double& cfl)
{
  double h,oneByAbsdt,density,viscosity,nrm_v;
  h = hFactor*elementDiameter;
  density = dm;
  viscosity =  a;
  nrm_v=0.0;
  for(int I=0;I<nSpace;I++)
    nrm_v+=f[I]*f[I];
  nrm_v = sqrt(nrm_v);
  cfl = nrm_v/h;
  oneByAbsdt =  fabs(dmt);
  tau0 = 1.0/(4.0*viscosity/(h*h) + 2.0*density*nrm_v/h + oneByAbsdt);
  tau1 = 4.0*viscosity + 2.0*density*nrm_v*h + oneByAbsdt*h*h;
}

inline
void calculateSubgridError_tau(int eN,
                               int k,
                               double  hFactor,
                               double* elementDiameter,
                               double* dmt,
                               double* dm,
                               double* f,
                               double* a,
                               double* tau0,
                               double* tau1,
                               double* cfl)
{
  calculateSubgridError_tau_c(hFactor,
			      elementDiameter[eN],
			      dmt[eN*nQuadraturePoints_element+k],
			      dm[eN*nQuadraturePoints_element+
				 k],
			      &f[eN*nQuadraturePoints_element*nSpace+
				k*nSpace],
			      a[eN*nQuadraturePoints_element*nSpace + 
				k*nSpace+1],
			      tau0[eN*nQuadraturePoints_element+k],
			      tau1[eN*nQuadraturePoints_element+k],
			      cfl[eN*nQuadraturePoints_element+k]);
}

inline
void calculateSubgridError_tauRes_c(const double& tau0,
				    const double& tau1,
				    const double& pdeResidualP,
				    const double& pdeResidualU,
				    const double& pdeResidualV,
				    const double& pdeResidualW,
				    double& subgridErrorP,
				    double& subgridErrorU,
				    double& subgridErrorV,
				    double& subgridErrorW)
{
  /* GLS momentum */
  subgridErrorU = tau0*pdeResidualU;
  subgridErrorV = tau0*pdeResidualV;
  subgridErrorW = tau0*pdeResidualW;
  /* GLS pressure */
  subgridErrorP = tau1*pdeResidualP;
}

inline
void calculateSubgridErrorDerivatives_tauRes_c(const double& tau0,
					       const double& tau1,
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
      dsubgridErrorP_du[j] = tau1*dpdeResidualP_du[j];
      dsubgridErrorP_dv[j] = tau1*dpdeResidualP_dv[j];
      dsubgridErrorP_dw[j] = tau1*dpdeResidualP_dw[j];
      /* GLS  momentum*/
      /* u */
      dsubgridErrorU_dp[j] = tau0*dpdeResidualU_dp[j];
      dsubgridErrorU_du[j] = tau0*dpdeResidualU_du[j];
      /* v */
      dsubgridErrorV_dp[j] = tau0*dpdeResidualV_dp[j];
      dsubgridErrorV_dv[j] = tau0*dpdeResidualV_dv[j];
      /* w */
      dsubgridErrorW_dp[j] = tau0*dpdeResidualW_dp[j];
      dsubgridErrorW_dw[j] = tau0*dpdeResidualW_dw[j];
    }
}

inline
void calculateSubgridError_tauRes(int eN,
                                  int k,
                                  double* tau0,
                                  double* tau1,
                                  double* pdeResidualP,
                                  double* dpdeResidualP_du,
                                  double* dpdeResidualP_dv,
                                  double* dpdeResidualP_dw,
                                  double* pdeResidualU,
                                  double* dpdeResidualU_dp,
                                  double* dpdeResidualU_du,
                                  double* pdeResidualV,
                                  double* dpdeResidualV_dp,
                                  double* dpdeResidualV_dv,
                                  double* pdeResidualW,
                                  double* dpdeResidualW_dp,
                                  double* dpdeResidualW_dw,
                                  double* subgridErrorP,
                                  double* dsubgridErrorP_du,
                                  double* dsubgridErrorP_dv,
                                  double* dsubgridErrorP_dw,
                                  double* subgridErrorU,
                                  double* dsubgridErrorU_dp,
                                  double* dsubgridErrorU_du,
                                  double* subgridErrorV,
                                  double* dsubgridErrorV_dp,
                                  double* dsubgridErrorV_dv,
                                  double* subgridErrorW,
                                  double* dsubgridErrorW_dp,
                                  double* dsubgridErrorW_dw)
{
  calculateSubgridError_tauRes_c(tau0[eN*nQuadraturePoints_element+k],
				 tau1[eN*nQuadraturePoints_element+k],
				 pdeResidualP[eN*nQuadraturePoints_element+k],
				 pdeResidualU[eN*nQuadraturePoints_element+k],
				 pdeResidualV[eN*nQuadraturePoints_element+k],
				 pdeResidualW[eN*nQuadraturePoints_element+k],
				 subgridErrorP[eN*nQuadraturePoints_element+k],
				 subgridErrorU[eN*nQuadraturePoints_element+k],
				 subgridErrorV[eN*nQuadraturePoints_element+k],
				 subgridErrorW[eN*nQuadraturePoints_element+k]);
  calculateSubgridErrorDerivatives_tauRes_c(tau0[eN*nQuadraturePoints_element+k],
				 tau1[eN*nQuadraturePoints_element+k],
				 &dpdeResidualP_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualP_dv[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualP_dw[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualU_dp[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualU_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualV_dp[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualV_dv[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualW_dp[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dpdeResidualW_dw[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorP_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorP_dv[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorP_dw[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorU_dp[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorU_du[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorV_dp[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorV_dv[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorW_dp[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element],
				 &dsubgridErrorW_dw[eN*nQuadraturePoints_element*nDOF_trial_element+k*nDOF_trial_element]);
}

inline
void calculateNumericalDiffusion_c(const double& shockCapturingDiffusion,
				   const double& elementDiameter,
				   const double& strong_residual,
				   const double grad_u[nSpace],
				   double& numDiff)
{
  double h,
    num,
    den,
    n_grad_u;
  h = elementDiameter;
  n_grad_u = 0.0;
  for (int I=0;I<nSpace;I++)
    n_grad_u += grad_u[I]*grad_u[I];
  num = shockCapturingDiffusion*0.5*h*fabs(strong_residual);
  den = sqrt(n_grad_u) + 1.0e-8;
  numDiff = num/den;
}

inline
void calculateNumericalDiffusion(int eN,
				 int k,
				 double shockCapturingDiffusion,
				 double* elementDiameter,
				 double* strong_residual,
				 double* grad_u,
				 double* numDiff)
{
  calculateNumericalDiffusion_c(shockCapturingDiffusion,
				elementDiameter[eN],
				strong_residual[eN*nQuadraturePoints_element+k],
				&grad_u[eN*nQuadraturePoints_element*nSpace+k*nSpace],
				numDiff[eN*nQuadraturePoints_element+k]);
}

inline
double SubgridError_c(const double& error,
		      const double& Lstar_w_dV)
{
  return -error*Lstar_w_dV;
}

inline
void updateSubgridError(int eN,
			int k,
			int i,
			double* error,
			double* Lstar_w_dV,
			double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] +=  SubgridError_c(error[eN*nQuadraturePoints_element + 
								   k],
							     Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
									k*nDOF_test_element + 
									i]);
}

inline
double SubgridErrorJacobian_c(const double& derror,
				  const double& Lstar_w_dV)
{
  return -derror*Lstar_w_dV;
}

inline
void updateSubgridErrorJacobian(int eN,
				int k,
				int j,
				int i,
				double* derror,
				double* Lstar_w_dV,
				double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element + 
			 j] += SubgridErrorJacobian_c(derror[eN*nQuadraturePoints_element*nDOF_trial_element+
							     k*nDOF_trial_element+
							     j],
						      Lstar_w_dV[eN*nQuadraturePoints_element*nDOF_test_element + 
								 k*nDOF_test_element + 
								 i]);
}

inline
double NumericalDiffusion_c(const double& numDiff,
			    const double grad_u[nSpace],
			    const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    tmp +=  numDiff*grad_u[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateNumericalDiffusion(int eN,
			      int k,
			      int i,
			      double* numDiff,
			      double* grad_u,
			      double* grad_w_dV,
			      double* weak_residual)
{
  weak_residual[eN*nDOF_test_element + i] += NumericalDiffusion_c(numDiff[eN*nQuadraturePoints_element + 
									  k],
								  &grad_u[eN*nQuadraturePoints_element*nSpace + 
									  k*nSpace],
								  &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
									     k*nDOF_test_element*nSpace + 
									     i*nSpace]);
}

inline
double NumericalDiffusionJacobian_c(const double& numDiff,
					const double grad_v[nSpace],
					const double grad_w_dV[nSpace])
{
  double tmp=0.0;
  for (int I=0;I<nSpace;I++)
    tmp += numDiff*grad_v[I]*grad_w_dV[I];
  return tmp;
}

inline
void updateNumericalDiffusionJacobian(int eN,
				      int k,
				      int j,
				      int i,
				      double* numDiff,
				      double* grad_v,
				      double* grad_w_dV,
				      double* jacobian_weak_residual)
{
  jacobian_weak_residual[eN*nDOF_test_X_trial_element + 
			 i*nDOF_trial_element + 
			 j] += NumericalDiffusionJacobian_c(numDiff[eN*nQuadraturePoints_element + 
								    k],
							    &grad_v[eN*nQuadraturePoints_element*nDOF_trial_element*nSpace + 
								    k*nDOF_trial_element*nSpace + 
								    j*nSpace],
							    &grad_w_dV[eN*nQuadraturePoints_element*nDOF_test_element*nSpace + 
								       k*nDOF_test_element*nSpace + 
								       i*nSpace]);
}

inline
void exteriorNumericalAdvectiveFlux_c(const int& isDOFBoundary_p,
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
void exteriorNumericalAdvectiveFluxDerivatives_c(const int& isDOFBoundary_p,
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
void exteriorNumericalAdvectiveFlux(int ebNE, 
				    int k, 
				    int *isDOFBoundary_p,
				    int *isDOFBoundary_u,
				    int *isDOFBoundary_v,
				    int *isDOFBoundary_w,
				    int *isFluxBoundary_p,
				    int *isFluxBoundary_u,
				    int *isFluxBoundary_v,
				    int *isFluxBoundary_w,
				    double* n,
				    double* bc_p,
				    double* bc_f_mass,
				    double* bc_f_umom,
				    double* bc_f_vmom,
				    double* bc_f_wmom,
				    double* bc_flux_mass,
				    double* bc_flux_umom,
				    double* bc_flux_vmom,
				    double* bc_flux_wmom,
				    double* p,
				    double* f_mass,
				    double* f_umom,
				    double* f_vmom,
				    double* f_wmom,
				    double* df_mass_du,
				    double* df_mass_dv,
				    double* df_mass_dw,
				    double* df_umom_dp,
				    double* df_umom_du,
				    double* df_umom_dv,
				    double* df_umom_dw,
				    double* df_vmom_dp,
				    double* df_vmom_du,
				    double* df_vmom_dv,
				    double* df_vmom_dw,
				    double* df_wmom_dp,
				    double* df_wmom_du,
				    double* df_wmom_dv,
				    double* df_wmom_dw,
				    double* flux_mass,
				    double* flux_umom,
				    double* flux_vmom,
				    double* flux_wmom,
				    double* dflux_mass_du,
				    double* dflux_mass_dv,
				    double* dflux_mass_dw,
				    double* dflux_umom_dp,
				    double* dflux_umom_du,
				    double* dflux_umom_dv,
				    double* dflux_umom_dw,
				    double* dflux_vmom_dp,
				    double* dflux_vmom_du,
				    double* dflux_vmom_dv,
				    double* dflux_vmom_dw,
				    double* dflux_wmom_dp,
				    double* dflux_wmom_du,
				    double* dflux_wmom_dv,
				    double* dflux_wmom_dw,
				    double* ebqe_velocity)
{
  exteriorNumericalAdvectiveFlux_c(isDOFBoundary_p[ebNE*nQuadraturePoints_elementBoundary+k],
				   isDOFBoundary_u[ebNE*nQuadraturePoints_elementBoundary+k],
				   isDOFBoundary_v[ebNE*nQuadraturePoints_elementBoundary+k],
				   isDOFBoundary_w[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_p[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_u[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_v[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_w[ebNE*nQuadraturePoints_elementBoundary+k],
				   &n[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   bc_p[ebNE*nQuadraturePoints_elementBoundary+k],
				   &bc_f_mass[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &bc_f_umom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &bc_f_vmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &bc_f_wmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   bc_flux_mass[ebNE*nQuadraturePoints_elementBoundary+k],
				   bc_flux_umom[ebNE*nQuadraturePoints_elementBoundary+k],
				   bc_flux_vmom[ebNE*nQuadraturePoints_elementBoundary+k],
				   bc_flux_wmom[ebNE*nQuadraturePoints_elementBoundary+k],
				   p[ebNE*nQuadraturePoints_elementBoundary+k],
				   &f_mass[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &f_umom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &f_vmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &f_wmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_mass_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_mass_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_mass_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_dp[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_dp[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_dp[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   flux_mass[ebNE*nQuadraturePoints_elementBoundary+k],
				   flux_umom[ebNE*nQuadraturePoints_elementBoundary+k],
				   flux_vmom[ebNE*nQuadraturePoints_elementBoundary+k],
				   flux_wmom[ebNE*nQuadraturePoints_elementBoundary+k],
				   &ebqe_velocity[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace]);
  exteriorNumericalAdvectiveFluxDerivatives_c(isDOFBoundary_p[ebNE*nQuadraturePoints_elementBoundary+k],
				   isDOFBoundary_u[ebNE*nQuadraturePoints_elementBoundary+k],
				   isDOFBoundary_v[ebNE*nQuadraturePoints_elementBoundary+k],
				   isDOFBoundary_w[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_p[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_u[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_v[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary_w[ebNE*nQuadraturePoints_elementBoundary+k],
				   &n[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   bc_p[ebNE*nQuadraturePoints_elementBoundary+k],
				   &bc_f_mass[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &bc_f_umom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &bc_f_vmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &bc_f_wmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   bc_flux_mass[ebNE*nQuadraturePoints_elementBoundary+k],
				   bc_flux_umom[ebNE*nQuadraturePoints_elementBoundary+k],
				   bc_flux_vmom[ebNE*nQuadraturePoints_elementBoundary+k],
				   bc_flux_wmom[ebNE*nQuadraturePoints_elementBoundary+k],
				   p[ebNE*nQuadraturePoints_elementBoundary+k],
				   &f_mass[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &f_umom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &f_vmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &f_wmom[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_mass_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_mass_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_mass_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_dp[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_umom_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_dp[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_vmom_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_dp[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_du[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_dv[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   &df_wmom_dw[ebNE*nQuadraturePoints_elementBoundary*nSpace+k*nSpace],
				   dflux_mass_du[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_mass_dv[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_mass_dw[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_umom_dp[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_umom_du[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_umom_dv[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_umom_dw[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_vmom_dp[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_vmom_du[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_vmom_dv[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_vmom_dw[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_wmom_dp[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_wmom_du[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_wmom_dv[ebNE*nQuadraturePoints_elementBoundary+k],
				   dflux_wmom_dw[ebNE*nQuadraturePoints_elementBoundary+k]);
}

inline
void exteriorNumericalDiffusiveFlux_c(const double& eps,
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
void exteriorNumericalDiffusiveFlux(double eps,
				    double phi,
				    int ebNE,
				    int k,
				    int* rowptr,
				    int* colind,
				    int* isDOFBoundary,
				    int* isFluxBoundary,
				    double* n,
				    double* bc_a,
				    double* bc_u,
				    double* bc_flux,
				    double* a,
				    double* grad_phi,
				    double* u,
				    double* penalty,
				    double* flux)
{
  const int nnz=rowptr[nSpace];
  exteriorNumericalDiffusiveFlux_c(eps,
				   phi,
				   rowptr,
				   colind,
				   isDOFBoundary[ebNE*nQuadraturePoints_elementBoundary+k],
				   isFluxBoundary[ebNE*nQuadraturePoints_elementBoundary+k],
				   &n[ebNE*nQuadraturePoints_elementBoundary*nSpace+
				      k*nSpace],
				   &bc_a[ebNE*nQuadraturePoints_elementBoundary*nnz+
					 k*nnz],
				   bc_u[ebNE*nQuadraturePoints_elementBoundary+
					k],
				   bc_flux[ebNE*nQuadraturePoints_elementBoundary+k],
				   &a[ebNE*nQuadraturePoints_elementBoundary*nnz+
				      k*nnz],
				   &grad_phi[ebNE*nQuadraturePoints_elementBoundary*nSpace+
					     k*nSpace],
				   u[ebNE*nQuadraturePoints_elementBoundary+
				     k],
				   penalty[ebNE*nQuadraturePoints_elementBoundary+
					   k],
				   flux[ebNE*nQuadraturePoints_elementBoundary+k]);
}

inline
double ExteriorElementBoundaryFlux_c(const double& flux,
				     const double& w_dS)
{
  return flux*w_dS;
}

inline
void updateExteriorElementBoundaryFlux(int ebNE,
				       int eN_global,
				       int k,
                                       int i,
                                       double* flux,
                                       double* w_dS,
                                       double* residual)
{
  residual[eN_global*nDOF_test_element+
	   i] +=  ExteriorElementBoundaryFlux_c(flux[ebNE*nQuadraturePoints_elementBoundary+
						     k],
						w_dS[ebNE*nQuadraturePoints_elementBoundary*nDOF_test_element+
						     k*nDOF_test_element+
						     i]);
}

inline
double ExteriorNumericalAdvectiveFluxJacobian_c(const double& dflux_left,
						const double& v)
{
  return dflux_left*v;
}

inline
void updateExteriorNumericalAdvectiveFluxJacobian(int ebNE,
						  int k,
						  int j,
						  double* dflux_left,
						  double* v,
						  double* fluxJacobian)
{
  fluxJacobian[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+k*nDOF_trial_element+j] 
    += 
    ExteriorNumericalAdvectiveFluxJacobian_c(dflux_left[ebNE*nQuadraturePoints_elementBoundary+k],
					     v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+k*nDOF_trial_element+j]);
}

inline
double ExteriorNumericalDiffusiveFluxJacobian_c(const double& eps,
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

inline
void updateExteriorNumericalDiffusiveFluxJacobian(double eps,
						  double phi,
						  int ebNE,
						  int k,
						  int j,
						  int* rowptr,
						  int* colind,
						  int* isDOFBoundary,
						  double* n,
						  double* a,
						  double* v,
						  double* grad_v,
						  double* penalty,
						  double* fluxJacobian)
{
  int nnz=rowptr[nSpace];
  fluxJacobian[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+
	       k*nDOF_trial_element+
	       j] += ExteriorNumericalDiffusiveFluxJacobian_c(eps,
							      phi,
							      rowptr,
							      colind,
							      isDOFBoundary[ebNE*nQuadraturePoints_elementBoundary+k],
							      &n[ebNE*nQuadraturePoints_elementBoundary*nSpace+
								 k*nSpace],
							      &a[ebNE*nQuadraturePoints_elementBoundary*nnz+
								 k*nnz],
							      v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element+
								k*nDOF_trial_element+
								j],
							      &grad_v[ebNE*nQuadraturePoints_elementBoundary*nDOF_trial_element*nSpace+
								      k*nDOF_trial_element*nSpace+
								      j*nSpace+
								      j],
							      penalty[ebNE*nQuadraturePoints_elementBoundary+
								      k]);
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
				double hFactor,
				double shockCapturingDiffusion,
				int* p_l2g, int* vel_l2g,
				double* elementDiameter,
				double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				// double* p_trial, double* vel_trial,
				// double* p_grad_trial, double* vel_grad_trial,
				// double* p_test_dV, double* vel_test_dV,
				// double* p_grad_test_dV, double* vel_grad_test_dV,
				// double* vel_Hess_trial,double* vel_Hess_test_dV,
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
				// double* p_trial_ext,
				// double* vel_trial_ext,
				// double* p_grad_trial_ext,
				// double* vel_grad_trial_ext,
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
				// double* p_test_dS_ext,
				// double* vel_test_dS_ext,
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
				double hFactor,
				double shockCapturingDiffusion,
				int* p_l2g, int* vel_l2g,
				double* elementDiameter,
				double* p_dof, double* u_dof, double* v_dof, double* w_dof,
				// double* p_trial, double* vel_trial,
				// double* p_grad_trial, double* vel_grad_trial,
				// double* p_test_dV, double* vel_test_dV,
				// double* p_grad_test_dV, double* vel_grad_test_dV,
				// double* vel_Hess_trial,double* vel_Hess_test_dV,
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
				// double* p_trial_ext,
				// double* vel_trial_ext,
				// double* p_grad_trial_ext,
				// double* vel_grad_trial_ext,
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
				// double* p_test_dS_ext,
				// double* vel_test_dS_ext,
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
  void calculateVelocityAverage_RANS2PV2(int nExteriorElementBoundaries_global,
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
